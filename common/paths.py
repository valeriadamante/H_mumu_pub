#!/usr/bin/env python3
"""Storage locations of the analysis, from config/paths.yaml.

    from common.paths import path
    skim = path("skim")                      # pathlib.Path

    python3 common/paths.py                  # every key
    python3 common/paths.py skim             # one value
    eval "$(python3 common/paths.py --shell)"  # export HMUMU_<KEY> (env.sh does it)

Precedence: HMUMU_<KEY> environment variable, then the file named by
HMUMU_PATHS, then config/paths.yaml.  "{key}" expands to another key and
"{repo}" to this checkout.
"""
from __future__ import annotations

import os
import shlex
import sys
from functools import lru_cache
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
DEFAULTS = REPO / "config" / "paths.yaml"


def env_name(key: str) -> str:
    return "HMUMU_" + key.upper()


@lru_cache(None)
def paths() -> dict[str, str]:
    raw = dict(yaml.safe_load(DEFAULTS.read_text()) or {})
    override = os.environ.get("HMUMU_PATHS")
    if override:
        raw.update(yaml.safe_load(Path(override).read_text()) or {})
    for key in list(raw):
        if env_name(key) in os.environ:
            raw[key] = os.environ[env_name(key)]

    resolved: dict[str, str] = {}

    def resolve(key: str, stack: tuple[str, ...] = ()) -> str:
        if key in resolved:
            return resolved[key]
        if key in stack:
            raise ValueError(f"config/paths.yaml: circular reference {' -> '.join(stack + (key,))}")
        if key not in raw:
            raise KeyError(f"config/paths.yaml has no key {key!r}")
        value = str(raw[key])
        while "{" in value:
            start = value.index("{")
            end = value.index("}", start)
            name = value[start + 1:end]
            replacement = str(REPO) if name == "repo" else resolve(name, stack + (key,))
            value = value[:start] + replacement + value[end + 1:]
        resolved[key] = value.rstrip("/") or "/"
        return resolved[key]

    return {key: resolve(key) for key in raw}


def path(key: str) -> Path:
    return Path(paths()[key])


def export_environment() -> None:
    """Set HMUMU_<KEY> for child processes (shell campaign configs read them)."""
    for key, value in paths().items():
        os.environ.setdefault(env_name(key), value)


def main(argv: list[str]) -> int:
    if argv[:1] == ["--shell"]:
        # Keys the user set by hand stay as they are; the ones exported here
        # are listed so that the next `source env.sh` can drop them first.
        exported = []
        for key, value in paths().items():
            if env_name(key) not in os.environ:
                print(f"export {env_name(key)}={shlex.quote(value)}")
                exported.append(env_name(key))
        print(f"export HMUMU_EXPORTED_KEYS={shlex.quote(' '.join(exported))}")
        return 0
    if argv:
        print(paths()[argv[0]])
        return 0
    for key, value in paths().items():
        print(f"{key:22s} {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
