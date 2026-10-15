"""Apply explicit, case-specific assistant decisions to measured evidence.

The TSVs are authored during individual inspection. No default decision is
inferred for missing rows, and an automatic draft can never mark itself reviewed.
"""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re

from ..io import file_ref, read_json, write_json

LABELS = {"N": "nonpeptide", "P": "peptide_or_conjugate", "H": "chemical_review_hold", "B": "unresolved_graph"}
BENIGN = {"MOL2_PARSE_FAILED", "HISTORICAL_CSV_PEPTIDE_LABEL_CHECKED", "BETA_LACTAM_CLASSIFICATION_REVIEW"}


def apply_reviews(root, require_complete=False):
    root = Path(root).resolve()
    selection = read_json(root / "selection_manifest.json")["selected_ids"]
    decisions = {}
    for path in sorted((root / "reviews").glob("*.tsv")):
        for line, row in enumerate(csv.reader(path.open(encoding="utf-8"), delimiter="\t"), 1):
            if len(row) != 3 or row[0] not in selection or row[1] not in LABELS or len(row[2]) < 8:
                raise ValueError(f"Invalid explicit decision {path}:{line}")
            if row[0] in decisions:
                raise ValueError(f"Duplicate individual review: {row[0]}")
            decisions[row[0]] = (row[1], row[2], file_ref(path), line)
    if require_complete and set(decisions) != set(selection):
        raise ValueError(f"Missing {len(set(selection) - set(decisions))} individual decisions")
    records, changes, counts = [], [], Counter()
    for ident in selection:
        path = root / "complexes" / ident / "curation_record.json"
        r = read_json(path)
        if ident not in decisions:
            counts["pending"] += 1
            continue
        code, note, ref, line = decisions[ident]
        if (code == "B") != (r["smiles"] is None):
            raise ValueError(f"{ident}: graph availability and review decision disagree")
        prior = r["peptide_analysis"]["status"]
        r["individual_review"] = {
            "status": "completed" if code != "B" else "blocked_no_graph",
            "decision": LABELS[code],
            "structure_reasoning_zh": note,
            "reviewer_kind": "assistant",
            "human_reviewed": False,
            "method": "individual canonical graph and structural/coordinate evidence review; not a claim of manual 3D visual inspection",
            "review_source": ref,
            "review_source_line": line,
            "evidence_at_review": r["evidence"],
            "corrections_zh": "以此条结构判读替换自动草稿的骨架/肽类别；几何距离不升级为氢键或活性证据。",
        }
        r["peptide_class"] = LABELS[code]
        r["training_exclusion_reasons"] = sorted(set(r["source_quality_flags"]) - BENIGN)
        if code != "N":
            r["training_exclusion_reasons"].append(LABELS[code])
        r["training_eligible"] = code == "N" and not r["training_exclusion_reasons"]
        r["training_label_kind"] = "observed_structure_caption" if r["training_eligible"] else None
        caption = note.split("；")[0].split("。")[0]
        # Audit notes may mention another case or an exclusion decision; those
        # are not part of the conditioning text. Never include PDB IDs/SMILES.
        if re.search(r"\b[0-9][a-z0-9]{3}\b|解析|回退|隔离|同图|同化学|复合物|冲突", caption):
            caption = "由提供的结构测得以下组成"
        r["structural_caption_zh"] = caption + "。" + r["automatic_draft"]["structural_caption"]
        write_json(path, r)
        (path.parent / "annotation.txt").write_text(
            f"{ident} | 第 {r['selection_rank']} 条 | {LABELS[code]}\n"
            f"逐条图结构判读：{note}\n结构摘要：{r['structural_caption_zh']}\n"
            f"源结构 SMILES：{r['smiles']}\n"
            f"断键初筛：{prior}；最大连通α氨基酸组：{r['peptide_analysis']['largest_connected_amino_acid_group']}\n"
            f"质量标记：{', '.join(r['source_quality_flags']) or '无'}\n"
            f"训练可用：{r['training_eligible']}；排除原因：{', '.join(r['training_exclusion_reasons']) or '无'}\n"
            "证据级别：观察结构的描述/重建，非优化配对；未证实氢键、溶剂暴露或亲和力改善。\n"
            "复核方式：逐条分子图与测量证据审阅；不是人工专家签名，也不宣称逐张进行3D目视检查。\n",
            encoding="utf-8",
        )
        counts[LABELS[code]] += 1
        counts["training_eligible"] += r["training_eligible"]
        records.append(r)
        if (
            (code == "N" and prior != "nonpeptide")
            or (code == "P" and prior not in {"peptide", "peptide_conjugate"})
            or code == "H"
        ):
            changes.append({"complex_id": ident, "automatic": prior, "reviewed": LABELS[code], "reason": note})
    with (root / "reviewed_records.jsonl").open("w", encoding="utf-8") as handle:
        for r in records:
            handle.write(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n")
    write_json(root / "individual_corrections.json", changes)
    report = {
        "selected": len(selection),
        "decisions_written": len(decisions),
        "counts": dict(counts),
        "all_decisions_present": len(decisions) == len(selection),
        "graph_review_completed": sum(r["smiles"] is not None for r in records),
        "graph_review_blocked": sum(r["smiles"] is None for r in records),
    }
    write_json(root / "review_summary.json", report)
    print(json.dumps(report, ensure_ascii=False))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    apply_reviews(args.root, args.require_complete)


if __name__ == "__main__":
    main()
