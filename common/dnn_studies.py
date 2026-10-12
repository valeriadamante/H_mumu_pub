"""Shared process composition for the DNN output studies.

Reads config/dnn_studies.yaml so that binning optimization, the weighted /
unweighted performance comparison and the sensitivity-vs-nbins scan all agree
on which processes enter signal and background.
"""

import os
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(
    os.environ.get("ANALYSIS_PATH", str(REPO))
) / "config" / "dnn_studies.yaml"

_CACHE = {}


def load_dnn_studies_config(path=None):
    """Return the parsed study configuration, failing loudly if absent."""
    config_path = Path(path) if path is not None else CONFIG_PATH
    key = str(config_path)
    if key in _CACHE:
        return _CACHE[key]
    if not config_path.is_file():
        raise FileNotFoundError(f"DNN study configuration not found: {config_path}")
    with config_path.open() as stream:
        config = yaml.safe_load(stream)
    for section in ("object", "signal", "background", "exclude", "binning"):
        if section not in config:
            raise KeyError(f"{config_path}: missing '{section}' section")
    _CACHE[key] = config
    return config


def signal_processes(generator="powheg", path=None):
    """Return the signal process names for one generator set."""
    signals = load_dnn_studies_config(path)["signal"]
    if generator == "all":
        return [name for key in sorted(signals) for name in signals[key]]
    if generator not in signals:
        raise KeyError(
            f"unknown signal set {generator!r}; available: {', '.join(sorted(signals))}"
        )
    return list(signals[generator])


def background_processes(path=None):
    return list(load_dnn_studies_config(path)["background"])


def exclude_patterns(path=None):
    return list(load_dnn_studies_config(path)["exclude"])


def object_path(region=None, variable=None, path=None):
    """Return '<region>/<variable>', defaulting to the configured pair."""
    configured = load_dnn_studies_config(path)["object"]
    return (
        f"{(region or configured['region']).rstrip('/')}"
        f"/{variable or configured['variable']}"
    )


def binning_constraints(path=None):
    return dict(load_dnn_studies_config(path)["binning"])
