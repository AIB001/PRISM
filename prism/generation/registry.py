"""Registry for the supported generation model identifiers."""

from typing import Any, Dict, Iterable, List, Mapping, Type

from .adapters import (
    DiffSBDDAdapter,
    FLOWRAdapter,
    MolCRAFTAdapter,
    Pocket2MolAdapter,
    PocketXMolAdapter,
    TargetDiffAdapter,
)
from .adapters.base import ModelAdapter
from .errors import ConfigurationError


ADAPTER_TYPES: Dict[str, Type[ModelAdapter]] = {
    "pocket2mol": Pocket2MolAdapter,
    "targetdiff": TargetDiffAdapter,
    "molcraft": MolCRAFTAdapter,
    "flowr": FLOWRAdapter,
    "pocketxmol": PocketXMolAdapter,
    "diffsbdd": DiffSBDDAdapter,
}


def model_ids() -> List[str]:
    return list(ADAPTER_TYPES)


def model_capabilities() -> Dict[str, Dict[str, Any]]:
    """Return stable, configuration-independent CLI capability metadata."""
    return {
        model: {
            "generation_mode": adapter_type.generation_mode,
            "requires_reference_ligand": adapter_type.requires_reference_ligand,
            "accepted_reference_ligand_suffixes": sorted(
                adapter_type.accepted_reference_ligand_suffixes
            ),
            "accepted_pocket_kinds": sorted(adapter_type.accepted_pocket_kinds),
        }
        for model, adapter_type in ADAPTER_TYPES.items()
    }


def expand_models(values: Iterable[str]) -> List[str]:
    requested = list(values)
    if "all" in requested:
        if len(requested) != 1:
            raise ConfigurationError("'all' cannot be combined with individual model names")
        return model_ids()
    unknown = sorted(set(requested) - set(ADAPTER_TYPES))
    if unknown:
        raise ConfigurationError(
            "Unknown generation model(s): " + ", ".join(unknown) + ". Available: " + ", ".join(model_ids())
        )
    return list(dict.fromkeys(requested))


def build_adapters(config: Mapping[str, Any], models: Iterable[str]) -> Dict[str, ModelAdapter]:
    model_configs = config.get("models", {})
    if not isinstance(model_configs, dict):
        raise ConfigurationError("Generation config 'models' must be a YAML mapping")
    result = {}
    for model in models:
        model_config = model_configs.get(model, {})
        if not isinstance(model_config, dict):
            raise ConfigurationError(f"Configuration for model '{model}' must be a YAML mapping")
        result[model] = ADAPTER_TYPES[model](model_config)
    return result
