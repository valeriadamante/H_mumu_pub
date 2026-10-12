#!/usr/bin/env python3
"""Normalize legacy scalar cross-section uncertainties into named mappings.

The transformation is intentionally line based so comments, entry order and
cross-section formulas remain untouched. It is idempotent: already-normalized
``unc`` mappings are preserved.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import yaml


LABELS = {
    "qcd scale": "qcd_scale",
    "scale": "scale",
    "theory": "theory",
    "th gaussian": "th_gaussian",
    "pdf+alpha s": "pdf_alpha_s",
    "pdf alpha s": "pdf_alpha_s",
    "pdf+alphas": "pdf_alpha_s",
    "pdf alphas": "pdf_alpha_s",
    "pdf+alphas": "pdf_alpha_s",
    "pdf alpha": "pdf_alpha_s",
    "pdf": "pdf",
    "alpha s": "alpha_s",
    "alphas": "alpha_s",
    "alpha": "alpha_s",
    "mq": "m_q",
    "mass": "mass",
    "ebeam": "beam_energy",
    "integration": "integration",
    "filter eff": "filter_efficiency",
    "total": "total",
    "e-nu": "electron_neutrino",
    "e+nu": "positron_neutrino",
}


def component_name(label: str) -> str:
    normalized = " ".join(label.lower().replace("_", " ").split())
    return LABELS.get(normalized, re.sub(r"[^a-z0-9]+", "_", normalized).strip("_"))


def yaml_scalar(raw: str) -> str:
    value: Any = yaml.safe_load(raw)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return raw
    return json.dumps(str(value), ensure_ascii=False)


def split_components(raw: str) -> list[tuple[str, str, bool]]:
    """Return (name, value, is_percentage) components from one legacy value."""
    groups = [group.strip() for group in raw.split(";")]
    components: list[tuple[str, str, bool]] = []
    used: set[str] = set()
    for group_index, group in enumerate(groups):
        suffix = "_BR" if group_index > 0 else ""
        cursor = 0
        found = False
        for match in re.finditer(r"\(([^()]*)\)", group):
            value = group[cursor:match.start()].strip(" ;")
            label = match.group(1).strip()
            cursor = match.end()
            # Parenthetical prose after a component is a note, not a new value.
            if not re.search(r"\d", value):
                continue
            value = re.sub(r"\s+for\s+(?:XS|BR)\s*$", "", value, flags=re.I).strip()
            name = component_name(label) + suffix
            if not name:
                continue
            base = name
            counter = 2
            while name in used:
                name = f"{base}_{counter}"
                counter += 1
            used.add(name)
            components.append((name, value.replace("%%", "%"), "%" in value))
            found = True
        if not found and re.search(r"\d", group):
            name = "total" + suffix
            components.append((name, group.replace("%%", "%"), "%" in group))
    if not components:
        return [("total", raw, "%" in raw)]
    return components


def uncertainty_lines(indent: str, raw: str, comment: str) -> list[str]:
    components = split_components(raw)
    lines = [f"{indent}unc:\n"]
    for index, (name, value, is_percentage) in enumerate(components):
        lines.append(f"{indent}  {name}:\n")
        suffix = f"  {comment}" if comment and index == 0 else ""
        lines.append(f"{indent}    value: {yaml_scalar(value)}{suffix}\n")
        if is_percentage:
            lines.append(f"{indent}    isPercentage: True\n")
    return lines


def split_comment(value: str) -> tuple[str, str]:
    # Existing uncertainty comments use a whitespace-delimited '#'.
    match = re.match(r"^(.*?)(?:\s{2,}(#.*))?$", value.rstrip())
    return match.group(1).strip(), (match.group(2) or "")


def normalize_text(text: str) -> str:
    lines = text.splitlines(keepends=True)
    output: list[str] = []
    in_entry = False
    entry_has_unc = False

    def finish_entry() -> None:
        nonlocal entry_has_unc
        if in_entry and not entry_has_unc:
            # Insert before trailing blank lines to keep each field in its block.
            trailing: list[str] = []
            while output and not output[-1].strip():
                trailing.append(output.pop())
            output.append("  unc: {}\n")
            output.extend(reversed(trailing))

    index = 0
    while index < len(lines):
        line = lines[index]
        if re.match(r"^[^\s#][^:]*:\s*(?:#.*)?$", line):
            finish_entry()
            in_entry = True
            entry_has_unc = False

        match = re.match(r"^(\s+)unc:\s*(.*?)\s*$", line)
        if not match:
            output.append(line)
            index += 1
            continue

        indent, rest = match.groups()
        entry_has_unc = True
        next_nonempty = index + 1
        while next_nonempty < len(lines) and not lines[next_nonempty].strip():
            next_nonempty += 1
        has_children = (
            next_nonempty < len(lines)
            and len(lines[next_nonempty]) - len(lines[next_nonempty].lstrip()) > len(indent)
        )

        if (not rest or rest in ("null", "{}")) and has_children:
            # Preserve an existing normalized mapping. The ``rest == '{}'``
            # case also repairs output produced by versions before the
            # idempotency check was added.
            output.append(f"{indent}unc:\n")
        elif not rest or rest == "null":
            output.append(f"{indent}unc: {{}}\n")
        elif rest.startswith("{") or rest.startswith("["):
            output.append(line)
        else:
            raw, comment = split_comment(rest)
            output.extend(uncertainty_lines(indent, raw, comment))
        index += 1

    finish_entry()
    return "".join(output)


def validate_schema(path: Path) -> None:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    for name, info in data.items():
        if not isinstance(info, dict) or not isinstance(info.get("unc"), dict):
            raise ValueError(f"{name}: unc is not a mapping")
        for component, detail in info["unc"].items():
            if not isinstance(detail, dict) or "value" not in detail:
                raise ValueError(f"{name}.unc.{component}: expected a mapping with value")
            if "isPercentage" in detail and detail["isPercentage"] is not True:
                raise ValueError(f"{name}.unc.{component}.isPercentage must be true")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", type=Path, default=Path("config/crossSections13p6TeV.yaml"))
    args = parser.parse_args()
    original = args.path.read_text(encoding="utf-8")
    normalized = normalize_text(original)
    args.path.write_text(normalized, encoding="utf-8")
    validate_schema(args.path)
    print(f"Normalized uncertainties in {args.path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
