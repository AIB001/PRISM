"""Bounded chat-completions client for DeepSeek, Qwen and local servers.

No provider SDK, model download, shell execution, or implicit network retry.
Only the explicitly configured environment variable supplies credentials.
"""

import ipaddress
import json
import os
import socket
from pathlib import Path
from urllib import error, parse, request

from .errors import DesignError, require
from .io import parse_json, read_json


def validate_json_schema(data, name):
    from jsonschema import Draft202012Validator

    schema = read_json(Path(__file__).parent / "schemas" / name)
    # These schemas only use local fragments. Older jsonschema resolvers lose
    # the document base for nested fragment references under a URN $id.
    schema.pop("$id", None)
    problems = sorted(Draft202012Validator(schema).iter_errors(data), key=lambda e: str(list(e.path)))
    if problems:
        # Do not echo values from model/config objects into error logs.
        fields = [".".join(map(str, p.path)) or "root" for p in problems[:5]]
        raise DesignError("Schema violation at: " + ", ".join(fields), "PLANNER_SCHEMA")


def load_llm_config(path):
    config = read_json(path)
    validate_json_schema(config, "planner-config-v1.json")
    url = parse.urlsplit(config["base_url"])
    require(
        bool(url.hostname) and not (url.username or url.password or url.query or url.fragment),
        "Planner base_url must not contain credentials, query parameters or fragments",
        "LLM_CONFIG",
    )
    try:
        local = url.hostname == "localhost" or ipaddress.ip_address(url.hostname).is_loopback
    except ValueError:
        local = url.hostname == "localhost"
    require(
        url.scheme == "https" or (local and url.scheme == "http"),
        "Use HTTPS, or HTTP on a loopback server",
        "LLM_CONFIG",
    )
    require(
        not url.path.rstrip("/").endswith("chat/completions"),
        "base_url should end before /chat/completions",
        "LLM_CONFIG",
    )
    require(
        local or config["api_key_env"] is not None,
        "A remote endpoint requires an explicit API-key environment variable",
        "LLM_CONFIG",
    )
    defaults = {
        "timeout_seconds": 60,
        "max_output_tokens": 4096,
        "max_attempts": 1,
        "temperature": 0.2,
        "json_mode": True,
        "request_options": {},
    }
    defaults.update(config)
    for field in ("max_output_tokens", "max_attempts"):
        require(type(defaults[field]) is int, field + " must be an integer", "LLM_CONFIG")
    require(
        not ({"thinking", "enable_thinking"} <= set(defaults["request_options"])),
        "Choose the thinking option supported by your provider",
        "LLM_CONFIG",
    )
    defaults["base_url"] = defaults["base_url"].rstrip("/")
    return defaults


class _NoRedirect(request.HTTPRedirectHandler):
    # Preserve the urllib callback signature even though every redirect is rejected.
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ARG002
        raise DesignError("Planner endpoint redirected; configure its final URL explicitly", "LLM_REDIRECT")


def chat_completion(config, messages):
    env_name = config["api_key_env"]
    key = os.environ.get(env_name) if env_name else None
    require(
        not env_name or bool(key),
        "Missing planner API key in configured environment variable",
        "LLM_CREDENTIALS_MISSING",
    )
    payload = {
        "model": config["model"],
        "messages": messages,
        "stream": False,
        "max_tokens": config["max_output_tokens"],
        "temperature": config["temperature"],
    }
    if config["json_mode"]:
        payload["response_format"] = {"type": "json_object"}
    payload.update(config["request_options"])
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if key:
        headers["Authorization"] = "Bearer " + key
    req = request.Request(
        config["base_url"] + "/chat/completions",
        data=json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8"),
        headers=headers,
        method="POST",
    )
    try:
        with request.build_opener(_NoRedirect()).open(req, timeout=config["timeout_seconds"]) as response:
            raw = response.read(2 * 1024 * 1024 + 1)
        require(len(raw) <= 2 * 1024 * 1024, "Planner HTTP response is too large", "LLM_RESPONSE_SIZE")
        data = parse_json(raw.decode("utf-8"))
    except error.HTTPError as exc:
        # Provider errors may echo request credentials or private prompts.
        raise DesignError(f"Planner HTTP error {exc.code}; no automatic network retry", "LLM_HTTP_ERROR") from None
    except (error.URLError, TimeoutError, socket.timeout, UnicodeError, OSError):
        raise DesignError(
            "Planner connection failed or timed out; no automatic network retry", "LLM_CONNECTION_ERROR"
        ) from None
    require(
        isinstance(data, dict) and isinstance(data.get("choices"), list) and len(data["choices"]) == 1,
        "Expected one chat-completion choice",
        "LLM_RESPONSE_FORMAT",
    )
    choice = data["choices"][0]
    require(
        isinstance(choice, dict) and isinstance(choice.get("message"), dict),
        "Invalid chat-completion message",
        "LLM_RESPONSE_FORMAT",
    )
    message = choice["message"]
    require(
        choice.get("finish_reason") == "stop" and not message.get("tool_calls") and not message.get("refusal"),
        "Planner response was incomplete, refused, or requested a tool call",
        "LLM_INCOMPLETE",
    )
    content = message.get("content")
    require(isinstance(content, str) and bool(content.strip()), "Planner returned no text", "LLM_EMPTY_RESPONSE")
    metadata = {
        "requested_model": config["model"],
        "returned_model": data.get("model"),
        "response_id": data.get("id"),
        "usage": data.get("usage"),
        "finish_reason": choice["finish_reason"],
        "provider_created": data.get("created"),
        "system_fingerprint": data.get("system_fingerprint"),
        "base_url": config["base_url"],
        "separate_model_api_invoked": True,
    }
    return content, metadata
