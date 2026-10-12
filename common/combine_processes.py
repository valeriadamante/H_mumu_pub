"""Combine processes: the single definition of dataset -> physics process -> card column.

Reads config/combine_processes.yaml.  Used by the hadd (one file per Combine
process, so merge-syst needs no change), by the campaign checks and by both card
builders (DNN and VBF selector), so that a dataset is added or moved in one place.
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config" / "combine_processes.yaml"


@lru_cache(None)
def config() -> dict:
    return yaml.safe_load(CONFIG.read_text())


def card_processes() -> dict:
    """Ordered {Combine process: entry} of the card columns."""
    return dict(config()["processes"])


def all_entries() -> dict:
    """Card columns plus the region-dependent processes (EWK), for the theory labels."""
    return {**card_processes(), **(config().get("region_processes") or {})}


def signal_processes() -> list[str]:
    return [name for name, entry in card_processes().items() if entry.get("signal")]


def background_processes() -> list[str]:
    return [name for name, entry in card_processes().items() if not entry.get("signal")]


def legacy() -> dict[str, list[str]]:
    """Merged process -> the Combine processes replacing it."""
    return dict(config().get("legacy") or {})


def split_processes() -> dict[str, dict]:
    """Combine processes that are a subset (``match``) of a process_names process."""
    return {name: entry for name, entry in card_processes().items() if entry.get("match")}


def members(entry) -> list[str]:
    entry = entry or {}
    return list(entry.get("datasets") or []) + list(entry.get("sub_processes") or [])


def combine_datasets(process_datasets: dict[str, list[str]]) -> dict[str, list[str]]:
    """{split Combine process: datasets} from an era's {process: datasets}.

    Only the processes that do not exist in process_names.yaml (the split ones);
    processes that are already a process_names entry are produced as before."""
    out = {}
    for name, entry in split_processes().items():
        if name in process_datasets:
            continue
        parent = process_datasets.get(entry["process"])
        if not parent:
            continue
        selected = [d for d in parent if re.search(entry["match"], d)]
        if selected:
            out[name] = selected
    return out


def check_partition(process_datasets: dict[str, list[str]]) -> list[str]:
    """Problems of the legacy splits: every dataset of a merged process must go to exactly one part."""
    problems = []
    split = combine_datasets(process_datasets)
    for merged, parts in legacy().items():
        datasets = process_datasets.get(merged) or []
        for dataset in datasets:
            owners = [p for p in parts if dataset in split.get(p, [])]
            if len(owners) != 1:
                problems.append(f"{merged}/{dataset}: in {owners or 'no'} split process")
    return problems


def shape_label_process(process: str) -> str:
    """Process whose systematics.yaml qcd_scale/pdf label a card process uses."""
    entry = all_entries().get(process) or {}
    return entry.get("shape_label", process)


def theory_label(process: str, dataset: str) -> str | None:
    """QCDscale/pdfAlphaS label of one dataset of a card process (None: no nuisance)."""
    if process in legacy():
        # A merged process standing in for its split ones: each dataset keeps the
        # label of the physics process it belongs to.
        part = next((p for p in legacy()[process]
                     if re.search(card_processes()[p]["match"], dataset)), None)
        return theory_label(part, dataset) if part else None
    theory = (all_entries().get(process) or {}).get("theory")
    if isinstance(theory, dict):
        return next((label for pattern, label in theory.items() if re.search(pattern, dataset)), None)
    return theory


def br_label(process: str) -> str | None:
    return (all_entries().get(process) or {}).get("br")
