#!/usr/bin/env python3
"""Simple lumi mask helper for raw golden JSON files."""

import json
import os
from pathlib import Path

__all__ = ["apply_lumi_filter", "load_lumi_map"]


def load_lumi_map(path):
    # Handle relative paths with ANALYSIS_PATH
    analysis_path = Path(os.environ.get("ANALYSIS_PATH", "."))
    path = Path(path)
    if not path.is_absolute():
        path = analysis_path / path
    if not path.exists():
        raise FileNotFoundError(f"Lumi JSON file not found: {path}")

    with path.open("r") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError("Lumi JSON must contain a dictionary of run -> lumi ranges")

    lumi_map = {}
    for run_str, lumi_ranges in data.items():
        if not isinstance(run_str, str):
            raise ValueError("Run keys must be strings")
        if not isinstance(lumi_ranges, list):
            raise ValueError(f"Lumi ranges for run {run_str} must be a list")

        run = int(run_str)
        ranges = []
        for lumi_range in lumi_ranges:
            if not isinstance(lumi_range, list) or len(lumi_range) != 2:
                raise ValueError(
                    f"Each lumi range for run {run_str} must be a two-element list"
                )
            start, end = int(lumi_range[0]), int(lumi_range[1])
            ranges.append((start, end))

        ranges.sort()
        lumi_map[run] = ranges

    return lumi_map


def _build_cpp_function(function_name, lumi_map):
    lines = [
        f"bool {function_name}(unsigned int run, unsigned int lumi) {{",
        "  switch (run) {",
    ]

    for run, ranges in sorted(lumi_map.items()):
        lines.append(f"    case {run}:")
        if not ranges:
            lines.append("      return false;")
        else:
            conds = [f"(lumi >= {start} && lumi <= {end})" for start, end in ranges]
            lines.append(f"      return {' || '.join(conds)};")
        lines.append("      break;")

    lines.extend([
        "    default:",
        "      return false;",
        "  }",
        "}",
    ])

    return "\n".join(lines)


def apply_lumi_filter(df, lumi_json_path, run_column="run", lumi_column="luminosityBlock", description="LumiFilter"):
    lumi_map = load_lumi_map(lumi_json_path)
    function_name = f"passesLumiMask_{abs(hash(str(lumi_json_path))) % 1000000}"
    cpp_code = _build_cpp_function(function_name, lumi_map)

    import ROOT
    ROOT.gInterpreter.Declare(cpp_code)

    filter_expr = f"{function_name}({run_column}, {lumi_column})"
    return df.Filter(filter_expr, description)