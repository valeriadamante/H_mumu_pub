#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path
import re
import sys
import time
import warnings

warnings.filterwarnings(
    "ignore",
    message="The value of the smallest subnormal.*type is zero.*",
    category=UserWarning,
    module=r"numpy\.core\.getlimits",
)

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join("/tmp", os.environ.get("USER", "user"), "matplotlib"),
)

import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
from matplotlib.backends.backend_pdf import PdfPages

# =========================================================
# Global style
# =========================================================

plt.style.use(hep.style.CMS)

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])

import common.utilities as utilities
from common.jet_component_splitting import pu_hard_component_style
from common.dy_component_patch import root_file_has_patch
from common.utilities import (
    initialize_root_runtime,
    RebinHisto,
    findBinEntry,
    findNewBins,
    getNewBins,
    is_valid_histogram,
)
from plotting_tools.plotting_functions import make_stacked_plot

initialize_root_runtime()


# =========================================================
# Sample helpers
# =========================================================

def normalize_sample_name(name):
    """
    Normalise the sample name.

    Examples:
      DY.root -> DY
      /path/to/TT.root -> TT
      DY -> DY
    """
    return os.path.splitext(os.path.basename(name))[0]


DY_012J_HARD_JETS = {
    "0J": 0,
    "1J PU": 0,
    "2J PU2": 0,
    "VBF PU2": 0,
    "1J Hard": 1,
    "2J PU1": 1,
    "VBF PU1": 1,
    "2J Hard": 2,
    "VBF Hard": 2,
}


def load_dy_012j_weights(era, weight_set="new"):
    era_name = str(era)
    if not era_name.startswith("Run3_"):
        era_name = f"Run3_{era_name}"
    if era_name.count("_") > 1:
        raise ValueError(
            "--dy-012j-reweight requires one physical era; a merged era such "
            f"as {era_name} has lost the per-era information"
        )
    repo = Path(os.environ["ANALYSIS_PATH"])
    if weight_set == "old":
        path = repo / "common" / "dy_component_patch_factors.json"
        payload = json.loads(path.read_text())
        era_text = era_name.removeprefix("Run3_")
        year = next((key for key in payload if era_text.startswith(key)), None)
        if year is None:
            raise ValueError(f"No old DY 0/1/2J weights for {era_name} in {path}")
        values = payload[year]
        hard0, hard1, hard2 = map(float, (values["2J PU"], values["1J PU"], values["2J Hard"]))
        weights = {
            "0J": hard0, "1JHard": hard1, "1JPU": hard0,
            "2JHard": hard2, "2JPU1": hard1, "2JPU2": hard0,
        }
        print(
            f"[INFO] Old DY 0/1/2J plot weights from {path}: "
            f"0J={hard0:.6g}, 1J={hard1:.6g}, 2J={hard2:.6g}"
        )
        return weights

    path = (
        repo / "reweights" / "dy_012j_reweight" / era_name
        / "dy_012j_reweight.json"
    )
    if not path.is_file():
        raise FileNotFoundError(f"DY 0/1/2J reweight JSON not found: {path}")
    payload = json.loads(path.read_text())
    corrections = payload.get("corrections", [])
    correction = next(
        (item for item in corrections if item.get("name") == "dy_012j_reweight"),
        None,
    )
    if correction is None:
        raise ValueError(f"Correction dy_012j_reweight not found in {path}")
    content = correction.get("data", {}).get("content", [])
    weights = {item["key"]: float(item["value"]) for item in content}
    expected = {
        "0J", "1JHard", "1JPU", "2JHard", "2JPU1", "2JPU2",
        "VBFHard", "VBFPU1", "VBFPU2",
    }
    if set(weights) != expected:
        raise ValueError(f"Expected six ggF and three VBF DY weights in {path}")
    print(f"[INFO] DY jet-component plot weights from {path}: {weights}")
    return weights


def dy_012j_plot_scale(sample_name, weights):
    if weights is None:
        return 1.0
    component = pu_hard_component_style(sample_name)
    if component is None or component[0] != "DY":
        return 1.0
    key = component[1].replace(" ", "")
    return weights.get(key, 1.0)


def parse_comma_separated_list(value):
    if value is None:
        return None

    items = []
    for item in value.split(","):
        item = item.strip()
        if item:
            items.append(item)

    return items


def classify_sample(sample_name, process_cfg):
    """
    Classify the sample using process_names.yaml.

    Rules:
      - is_data: True   -> data
      - is_signal: True -> signal
      - otherwise       -> background
    """
    if sample_name not in process_cfg:
        return None

    cfg = process_cfg[sample_name]

    if cfg.get("skip_plotting", False):
        return None

    is_data = bool(cfg.get("is_data", False))
    is_signal = bool(cfg.get("is_signal", False))

    if is_data and is_signal:
        raise RuntimeError(
            f"Sample {sample_name} is marked as both data and signal."
        )

    if is_data:
        sample_type = "data"
    elif is_signal:
        sample_type = "signal"
    else:
        sample_type = "background"

    return {
        "type": sample_type,
        "is_data": is_data,
        "is_signal": is_signal,
        "color": cfg.get(
            "color_mplhep",
            cfg.get("color", "black"),
        ),
        "name": cfg.get("name", sample_name),
    }


def expand_requested_samples(requested_samples, plot_groups_cfg):
    if requested_samples is None:
        return None

    expanded = set()
    background_groups = plot_groups_cfg.get("background_groups", {})
    other_group_key = plot_groups_cfg.get("other_group", {}).get("key", "OTHER")

    for sample in requested_samples:
        if sample == other_group_key:
            return None

        if sample in background_groups:
            expanded.update(get_group_members(background_groups[sample]))
        else:
            expanded.add(sample)

    return expanded


def get_group_members(group_cfg):
    members = []

    for key in ("processes", "sub_processes", "datasets", "aliases"):
        members.extend(group_cfg.get(key, []))

    return members


def get_group_color(group_cfg, fallback):
    return group_cfg.get("color_mplhep", group_cfg.get("color", fallback))


def get_group_member_info(plot_groups_cfg):
    member_info = {}

    for group_name, group_cfg in plot_groups_cfg.get("background_groups", {}).items():
        for member_name in get_group_members(group_cfg):
            member_info.setdefault(
                member_name,
                {
                    "group": group_name,
                    "group_cfg": group_cfg,
                },
            )

    return member_info


def get_process_scale_factors(plot_groups_cfg):
    scale_factors = {}

    for group_cfg in plot_groups_cfg.get("background_groups", {}).values():
        if "scale_factor" not in group_cfg:
            continue

        for member_name in get_group_members(group_cfg):
            scale_factors[member_name] = float(group_cfg["scale_factor"])

    for process_name, scale_cfg in plot_groups_cfg.get("process_scales", {}).items():
        if isinstance(scale_cfg, dict):
            scale_factors[process_name] = float(scale_cfg.get("scale_factor", 1.0))
        else:
            scale_factors[process_name] = float(scale_cfg)

    return scale_factors


def classify_plot_sample(sample_name, process_cfg, plot_groups_cfg):
    configured_component_styles = {
        process_name: process_info["jet_component_styles"]
        for process_name, process_info in process_cfg.items()
        if isinstance(process_info, dict)
        and process_info.get("jet_component_styles")
    }
    component_style = pu_hard_component_style(
        sample_name,
        configured_component_styles,
    )
    if component_style is not None:
        family, component_label, color = component_style
        if color is None:
            raise KeyError(
                f"Missing color for {family} {component_label} in "
                "process_names.yaml:jet_component_styles"
            )
        return {
            "type": "background",
            "is_data": False,
            "is_signal": False,
            "color": color,
            "name": f"{family} {component_label}",
        }

    sample_info = classify_sample(sample_name, process_cfg)

    if sample_info is not None and (
        sample_info["is_data"] or sample_info["is_signal"]
    ):
        return sample_info

    signal_styles = plot_groups_cfg.get("signal_styles", {})

    if sample_name in signal_styles:
        style_cfg = signal_styles[sample_name]
        return {
            "type": "signal",
            "is_data": False,
            "is_signal": True,
            "color": style_cfg.get("color_mplhep", style_cfg.get("color", "black")),
            "name": style_cfg.get("name", sample_name),
        }

    group_member_info = get_group_member_info(plot_groups_cfg)

    if sample_name in group_member_info:
        group_cfg = group_member_info[sample_name]["group_cfg"]
        return {
            "type": "background",
            "is_data": False,
            "is_signal": False,
            "color": get_group_color(group_cfg, "black"),
            "name": group_cfg.get("name", sample_name),
        }

    if sample_info is not None:
        return sample_info

    return None


def make_custom_sample_info(sample_name, index=0):
    colors = [
        "black",
        "red",
        "dodgerblue",
        "darkorange",
        "forestgreen",
        "purple",
    ]
    return {
        "type": "signal",
        "is_data": False,
        "is_signal": True,
        "color": colors[index % len(colors)],
        "name": sample_name,
    }


def make_group_process(group_name, group_cfg, members, input_processes):
    output_process = {
        "input": ",".join(input_processes[name]["input"] for name in members),
        "color": get_group_color(group_cfg, input_processes[members[0]]["color"]),
        "name": group_cfg.get("name", group_name),
        "is_data": False,
        "is_signal": False,
        "type": "background",
        "aliases": list(members),
        "systematic_inputs": [
            path
            for name in members
            for path in input_processes[name].get("systematic_inputs", [])
        ],
        "hists": {},
    }

    categories = {
        category
        for name in members
        for category in input_processes[name]["hists"]
    }
    for category in categories:
        member_hists = {
            name: input_processes[name]["hists"].get(category, {})
            for name in members
        }
        hist_names = {hist_name for hists in member_hists.values() for hist_name in hists}
        output_process["hists"].setdefault(category, {})

        for hist_name in hist_names:
            # For a shifted group template, members without that nuisance must
            # contribute their nominal histogram.  Otherwise the variation is
            # incorrectly normalized to only the affected subset of the group.
            nominal_name = re.sub(
                r"_(?:CMS_|QCD_|pdf_|JEReta).*(?:Up|Down)$", "", hist_name
            )
            is_shift = nominal_name != hist_name
            group_hist = None
            for member_name in members:
                hists = member_hists[member_name]
                hist = hists.get(hist_name)
                if hist is None and is_shift:
                    hist = hists.get(nominal_name)
                if hist is None:
                    continue
                if group_hist is None:
                    clone_name = f"{group_name}_{category}_{hist_name}"
                    group_hist = clone_hist_for_group(hist, clone_name)
                else:
                    group_hist.Add(hist)
            if group_hist is not None:
                output_process["hists"][category][hist_name] = group_hist

    return output_process


def clone_hist_for_group(hist, name):
    safe_name = name.replace("/", "_").replace(" ", "_")
    cloned = hist.Clone(safe_name)
    cloned.SetDirectory(0)
    return cloned


def apply_signal_styles(input_processes, plot_groups_cfg):
    signal_styles = plot_groups_cfg.get("signal_styles", {})

    for process_name, style_cfg in signal_styles.items():
        if process_name not in input_processes:
            continue

        if "name" in style_cfg:
            input_processes[process_name]["name"] = style_cfg["name"]

        if "color" in style_cfg:
            input_processes[process_name]["color"] = style_cfg["color"]
        if "color_mplhep" in style_cfg:
            input_processes[process_name]["color"] = style_cfg["color_mplhep"]


def apply_background_groups(input_processes, plot_groups_cfg, active_group_names=None):
    """
    Merge background macro-samples into plotting groups.

    Processes not listed in config/plot/process_groups.yaml are kept unchanged.
    Data and signals are never merged here.
    """
    background_groups = plot_groups_cfg.get("background_groups", {})

    if active_group_names is not None:
        background_groups = {
            group_name: group_cfg
            for group_name, group_cfg in background_groups.items()
            if group_name in active_group_names
        }

    if not background_groups:
        return input_processes

    output_processes = {}
    grouped_processes = set()

    for group_name, group_cfg in background_groups.items():
        member_names = get_group_members(group_cfg)
        members = [
            name
            for name in member_names
            if name in input_processes
            and name not in grouped_processes
            and not input_processes[name].get("is_data", False)
            and not input_processes[name].get("is_signal", False)
        ]

        if len(members) == 0:
            continue

        grouped_processes.update(members)

        output_processes[group_name] = make_group_process(
            group_name,
            group_cfg,
            members,
            input_processes,
        )

    other_group_cfg = plot_groups_cfg.get("other_group", {})
    other_name = other_group_cfg.get("key", "OTHER")
    use_other_group = (
        (
            active_group_names is None
            or other_name in active_group_names
        )
        and other_group_cfg.get("enabled", True)
    )

    if use_other_group:
        other_members = [
            process_name
            for process_name, process_info in input_processes.items()
            if process_name not in grouped_processes
            # PU/hard components must remain separate in the main stack.
            and pu_hard_component_style(process_name) is None
            and not process_info.get("is_data", False)
            and not process_info.get("is_signal", False)
        ]

        if len(other_members) > 0:
            grouped_processes.update(other_members)
            output_processes[other_name] = make_group_process(
                other_name,
                other_group_cfg,
                other_members,
                input_processes,
            )

    for process_name, process_info in input_processes.items():
        if process_name in grouped_processes:
            continue

        output_processes[process_name] = process_info

    return output_processes


def apply_plot_groups(input_processes, plot_groups_cfg, active_group_names=None):
    apply_signal_styles(input_processes, plot_groups_cfg)
    return apply_background_groups(
        input_processes,
        plot_groups_cfg,
        active_group_names=active_group_names,
    )


# =========================================================
# Helpers for histogram reading
# =========================================================

SYSTEMATIC_MARKERS = ("_CMS_", "_QCD_", "_pdf_", "_JEReta")


def nominal_of(hist_name):
    """The nominal histogram a shifted template belongs to (itself if nominal)."""
    if not hist_name.endswith(("Up", "Down")):
        return hist_name
    for marker in SYSTEMATIC_MARKERS:
        if marker in hist_name:
            return hist_name.split(marker, 1)[0]
    return hist_name


def sum_era_inputs(per_era, region_path):
    """Sum the processes read from several eras, without a merged file.

    Nominal: the sum over the eras.  Shifted template of a nuisance: the summed
    nominal plus, for each era that has it, (shifted - nominal) of that era.
    Nuisances are era-specific (CMS_pileup_2024, ...), so this is what the
    era-template merger writes after padding with the nominal, at the cost of a
    few histogram additions instead of rewriting every key in a file.
    """
    processes, found = {}, set()
    for era_processes, era_found in per_era:
        found |= era_found
        for name, info in era_processes.items():
            if name not in processes:
                processes[name] = dict(info, hists={region_path: {}}, systematic_inputs=[])
            processes[name]["systematic_inputs"] += info.get("systematic_inputs", [])
    for name, info in processes.items():
        era_hists = [p[name]["hists"].get(region_path, {}) for p, _ in per_era if name in p]
        out = info["hists"][region_path]
        for hists in era_hists:
            for key, hist in hists.items():
                if nominal_of(key) != key:
                    continue
                if key in out:
                    out[key].Add(hist)
                else:
                    out[key] = hist.Clone(f"{hist.GetName()}_eras")
                    out[key].SetDirectory(0)
        for hists in era_hists:
            for key, hist in hists.items():
                base = nominal_of(key)
                if base == key or base not in out:
                    continue
                if key not in out:
                    out[key] = out[base].Clone(f"{hist.GetName()}_eras")
                    out[key].SetDirectory(0)
                out[key].Add(hist)
                if base in hists:
                    out[key].Add(hists[base], -1.0)
    return processes, found


def get_available_histograms(
    root_file,
    region_path,
    scale_factor=1.0,
    recursive=True,
    exclude_2d=True,
):
    """
    Returns:
        [(histogram, hist_name), ...]
    """
    output = []

    directory = root_file.Get(region_path)

    if not directory:
        return output

    def scan_dir(tdir, prefix=""):
        for key in tdir.GetListOfKeys():
            name = key.GetName()
            object_path = f"{region_path}/{prefix}{name}"

            try:
                obj = key.ReadObj()
            except Exception as exc:
                print(
                    f"[WARNING] Cannot read {object_path} from "
                    f"{root_file.GetName()}: {exc}. Skip."
                )
                continue

            # A damaged ROOT key can make ReadObj return a null PyROOT proxy
            # instead of raising (often together with an R__unzip_header error).
            if not obj:
                print(
                    f"[WARNING] Null or corrupted ROOT object: {object_path} "
                    f"in {root_file.GetName()}. Skip."
                )
                continue

            # TH2/TH3 inherit from TH1 in ROOT, so InheritsFrom("TH1") alone
            # is not sufficient to select objects supported by this 1D plotter.
            if obj.InheritsFrom("TH1"):
                if exclude_2d and obj.GetDimension() != 1:
                    continue

                hist_name = f"{prefix}{name}" if prefix else name

                hist = obj.Clone(f"{hist_name}_clone")
                hist.SetDirectory(0)
                hist.Scale(scale_factor)

                output.append((hist, hist_name))

            elif recursive and obj.InheritsFrom("TDirectory"):
                new_prefix = f"{prefix}{name}/" if prefix else f"{name}/"
                scan_dir(obj, new_prefix)

    scan_dir(directory)

    return output


# =========================================================
# Main
# =========================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--era",
        required=True,
        type=str,
        help="Era, e.g. 2022, 2022EE, 2023, 2023BPix",
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Directory containing the ROOT files of the macro-samples",
    )

    parser.add_argument(
        "--output",
        default="plots_output",
        type=str,
        help="Output directory for plots",
    )

    parser.add_argument(
        "--region",
        default="Z_sideband_baseline",
        type=str,
        help="Region to plot",
    )

    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help=(
            "List of macro-samples to process, e.g.: "
            "DY TT SingleTop Data ggH. "
            "The name must match the ROOT file without .root, "
            "a key in process_names.yaml, or a group "
            "in config/plot/process_groups.yaml."
        ),
    )
    parser.add_argument(
        "--sample-color",
        action="append",
        default=[],
        metavar="SAMPLE=COLOR",
        help=(
            "Override the plotting color of one loaded sample/group. Repeat "
            "for multiple samples, for example --sample-color Flash=red."
        ),
    )
    parser.add_argument(
        "--sample-label",
        action="append",
        default=[],
        metavar="SAMPLE=LABEL",
        help=(
            "Override the legend label of one loaded sample/group. Repeat "
            "for multiple samples."
        ),
    )
    parser.add_argument(
        "--sample-input",
        action="append",
        default=[],
        metavar="SAMPLE=PATH",
        help=(
            "Load one sample from an alternate ROOT file or directory. "
            "Repeat for multiple samples; directory values resolve to "
            "PATH/SAMPLE.root."
        ),
    )

    parser.add_argument(
        "--systematics",
        action="store_true",
        help="Include systematic uncertainties",
    )
    parser.add_argument(
        "--systematic-input",
        action="append",
        default=[],
        metavar="NAME=DIR",
        help=(
            "Load shifted templates from a separate hadded directory while "
            "keeping nominal templates in --input; repeat for each family."
        ),
    )

    parser.add_argument(
        "--systematicGroup",
        action="append",
        default=None,
        help="Only draw this systematic group; repeat for multiple groups",
    )

    parser.add_argument(
        "--overlaySystematic",
        action="store_true",
        help="Overlay nominal and systematic Up/Down in the main panel",
    )

    parser.add_argument(
        "--totalSystematics",
        action="store_true",
        help="Draw the quadrature sum of all displayed systematic groups",
    )

    parser.add_argument(
        "--noMCStatUncertainty",
        action="store_true",
        help="Hide the MC statistical uncertainty band in the ratio panel",
    )

    parser.add_argument(
        "--wantData",
        action="store_true",
        help="Include data in plots and draw ratio",
    )

    parser.add_argument(
        "--wantLogY",
        action="store_true",
        help="Set y-axis to log scale",
    )

    parser.add_argument(
        "--logUncertainties",
        action="store_true",
        help="Use a logarithmic y-axis in the ratio/uncertainty panel",
    )

    parser.add_argument(
        "--rebin",
        action="store_true",
        help="Rebin histograms",
    )

    stack_group = parser.add_mutually_exclusive_group()
    stack_group.add_argument(
        "--stack",
        dest="do_stack",
        action="store_true",
        help="Draw backgrounds in a stack (default)",
    )
    stack_group.add_argument(
        "--no-stack",
        dest="do_stack",
        action="store_false",
        help="Draw backgrounds overlaid instead of stacked",
    )
    parser.set_defaults(do_stack=True)

    fill_group = parser.add_mutually_exclusive_group()
    fill_group.add_argument(
        "--fill-hists",
        dest="fill_hists",
        action="store_true",
        help="Draw non-line histograms filled (default)",
    )
    fill_group.add_argument(
        "--no-fill-hists",
        dest="fill_hists",
        action="store_false",
        help="Draw backgrounds as unfilled step lines",
    )
    parser.set_defaults(fill_hists=True)

    parser.add_argument(
        "--ratio-reference",
        default=None,
        type=str,
        help=(
            "Sample/group to use as ratio denominator. "
            "Default is Data/MC when --wantData is used. "
            "Example: --ratio-reference EWK_Herwig"
        ),
    )
    parser.add_argument(
        "--allow-custom-samples",
        action="store_true",
        help=(
            "Allow ROOT file names not present in process_names.yaml or "
            "process_groups.yaml. Useful for direct file-to-file comparisons."
        ),
    )

    parser.add_argument(
        "--normalize-dy-to-data",
        action="store_true",
        help=(
            "Scale one DY sample/group in each plot so that its integral "
            "matches the data integral."
        ),
    )

    parser.add_argument(
        "--normalize-mc-to-data",
        action="store_true",
        help=(
            "Scale all background MC samples in each plot with one common "
            "factor so that the sum of MC integrals matches the data integral."
        ),
    )

    parser.add_argument(
        "--dy-normalization-sample",
        default="DY",
        type=str,
        help=(
            "Sample(s) or plotting group(s) to scale with one common factor "
            "when --normalize-dy-to-data is used. Separate multiple targets "
            "with commas. Default: DY; if no literal DY background is present, "
            "all background samples/groups whose key or label starts with DY "
            "are scaled. Examples: DY, DY_amcatnlo, "
            "DYto2Mu_MLL105To160_ptll,DYto2Mu_MLL105To160_VBFFiltered_ptll."
        ),
    )
    parser.add_argument(
        "--component-composition",
        "--dy-composition",
        "--dy-ewk-composition",
        dest="dy_composition",
        action="store_true",
        help=(
            "Automatically add one per-bin component-fraction panel for every "
            "process family represented by at least two component samples. "
            "The DY/EWK option names remain as backward-compatible aliases."
        ),
    )

    parser.add_argument(
        "--dy-012j-weights",
        choices=("none", "old", "new"),
        default="none",
        help=(
            "DY component weights applied only while plotting: none (default), "
            "old (common/dy_component_patch_factors.json), or new (era correction JSON)"
        ),
    )
    parser.add_argument(
        "--dy-012j-reweight",
        dest="dy_012j_weights",
        action="store_const",
        const="new",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--vars",
        "--variables",
        default=None,
        type=str,
        help=(
            "Comma-separated list of variables to plot, "
            "for example: --vars m_mumu,DNN_NNOutput. "
            "If omitted, all variables found in the selected region are plotted."
        ),
    )

    parser.add_argument(
        "--exclude-2d",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Skip TH2/TH3 objects (default). Use --no-exclude-2d only for "
            "specialized workflows that handle multidimensional histograms."
        ),
    )
    parser.add_argument(
        "--multipage-pdf",
        nargs="?",
        const="all_plots.pdf",
        default=None,
        metavar="FILE",
        help=(
            "Save every variable as one page of a single PDF instead of "
            "writing separate PNG and PDF files. If FILE is omitted, use "
            "all_plots.pdf inside the era/region output directory."
        ),
    )


    parser.add_argument("--combined-eras", help="Comma-separated physical eras contributing to this merged plot")
    parser.add_argument("--x-rebin", help="Comma-separated bin edges for the plotted variables, "
                        "overriding x_rebin of config/plot/histograms.yaml (needs --rebin)")
    parser.add_argument("--era-inputs", help="Comma-separated per-era input directories, summed at plot time "
                        "(nominal summed; each era-specific variation = summed nominal + its shift); "
                        "replaces --input; pass the same eras to --combined-eras")
    args = parser.parse_args()
    combined_eras = [e if e.startswith("Run3_") else "Run3_" + e for e in (args.combined_eras or "").split(",") if e]
    if combined_eras and len(set(combined_eras)) != len(combined_eras):
        parser.error("--combined-eras contains duplicates")

    def parse_sample_overrides(values, option_name):
        result = {}
        for value in values:
            if "=" not in value:
                parser.error(f"{option_name} expects SAMPLE=VALUE, got: {value}")
            sample, override = value.split("=", 1)
            if not sample or not override:
                parser.error(f"{option_name} expects SAMPLE=VALUE, got: {value}")
            result[sample] = override
        return result

    sample_color_overrides = parse_sample_overrides(
        args.sample_color, "--sample-color"
    )
    sample_label_overrides = parse_sample_overrides(
        args.sample_label, "--sample-label"
    )
    sample_input_overrides = {}
    for sample, path in parse_sample_overrides(
        args.sample_input, "--sample-input"
    ).items():
        sample = normalize_sample_name(sample)
        if not path.endswith(".root"):
            path = os.path.join(path, f"{sample}.root")
        sample_input_overrides[sample] = os.path.abspath(path)
    if args.normalize_dy_to_data and args.normalize_mc_to_data:
        parser.error(
            "--normalize-dy-to-data and --normalize-mc-to-data are mutually exclusive"
        )

    try:
        dy_012j_weights = (
            load_dy_012j_weights(args.era, args.dy_012j_weights)
            if args.dy_012j_weights != "none" else None
        )
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))

    requested_variables = parse_comma_separated_list(args.vars)
    requested_variables_set = (
        set(requested_variables) if requested_variables is not None else None
    )

    startTime = time.time()

    # =====================================================
    # Configs
    # =====================================================

    cfg_dir = os.path.join(
        os.environ["ANALYSIS_PATH"],
        "config",
        args.era,
    )
    if combined_eras:
        cfg_dir = os.path.join(os.environ["ANALYSIS_PATH"], "config", combined_eras[-1])
    if not os.path.isdir(cfg_dir):
        combined_config_fallbacks = {
            "Run3_2022_25": "Run3_2025",
        }
        fallback_era = combined_config_fallbacks.get(args.era)
        if fallback_era is None:
            raise FileNotFoundError(
                f"No plotting configuration found for era {args.era}: "
                f"{cfg_dir}"
            )
        cfg_dir = os.path.join(
            os.environ["ANALYSIS_PATH"],
            "config",
            fallback_era,
        )
        print(
            f"[INFO] Using {fallback_era} process/selection definitions "
            f"for combined era {args.era}"
        )

    main_cfg = utilities.get_config(
        os.path.join(cfg_dir, "maincfg.yaml")
    )

    process_cfg = utilities.get_config(
        os.path.join(cfg_dir, "process_names.yaml")
    )

    plot_groups_cfg = utilities.get_config(
        os.path.join(
            os.environ["ANALYSIS_PATH"],
            "config",
            "plot",
            "process_groups.yaml",
        )
    )

    sel_cfg = utilities.get_config(
        os.path.join(cfg_dir, "selections.yaml")
    )

    syst_cfg = utilities.get_config(
        os.path.join(cfg_dir, "systematics.yaml")
    )

    hist_cfg = utilities.get_config(
        os.path.join(
            os.environ["ANALYSIS_PATH"],
            "config",
            "plot",
            "histograms.yaml",
        )
    )

    if args.x_rebin:
        edges = [float(edge) for edge in args.x_rebin.split(",")]
        for variable in (requested_variables or []):
            hist_cfg.setdefault(variable, {})["x_rebin"] = edges

    additional_cfg = utilities.get_config(
        os.path.join(
            os.environ["ANALYSIS_PATH"],
            "config",
            "plot",
            f"{combined_eras[-1] if combined_eras else args.era}.yaml",
        )
    )

    if combined_eras:
        luminosity = sum(float(utilities.get_config(os.path.join(os.environ["ANALYSIS_PATH"], "config", e, "maincfg.yaml"))["luminosity"]) for e in combined_eras)
        additional_cfg.setdefault("lumi_text", {})["text"] = f"{luminosity / 1000.0:.1f}"

    page_cfg = utilities.get_config(
        os.path.join(
            os.environ["ANALYSIS_PATH"],
            "config",
            "plot",
            "cms_stacked.yaml",
        )
    )

    region_path = args.region

    config_setup = {
        **page_cfg,
        **additional_cfg,
        **hist_cfg,
    }

    config_setup["wantLogY"] = args.wantLogY
    process_scale_factors = get_process_scale_factors(plot_groups_cfg)

    # =====================================================
    # Requested samples
    # =====================================================

    requested_samples = None
    requested_plot_groups = None

    if args.samples is not None:

        requested_samples_raw = set(
            normalize_sample_name(s)
            for s in args.samples
        )
        group_member_info = get_group_member_info(plot_groups_cfg)

        requested_plot_groups = {
            sample
            for sample in requested_samples_raw
            if (
                sample in plot_groups_cfg.get("background_groups", {})
                or sample == plot_groups_cfg.get("other_group", {}).get("key", "OTHER")
            )
        }
        requested_plot_groups.update(
            group_member_info[sample]["group"]
            for sample in requested_samples_raw
            if sample in group_member_info
        )
        if len(requested_plot_groups) == 0:
            requested_plot_groups = None

        requested_samples = expand_requested_samples(
            requested_samples_raw,
            plot_groups_cfg,
        )

        print("\nRequested macro-samples:")


        for sample in sorted(requested_samples_raw):

            if (
                sample not in process_cfg
                and sample not in plot_groups_cfg.get("background_groups", {})
                and sample != plot_groups_cfg.get("other_group", {}).get("key", "OTHER")
                and sample not in get_group_member_info(plot_groups_cfg)
                and sample not in plot_groups_cfg.get("signal_styles", {})
            ):
                component_info = classify_plot_sample(
                    sample, process_cfg, plot_groups_cfg
                )
                if component_info is not None:
                    print(f"  {sample}: {component_info['type']}")
                    continue
                if args.allow_custom_samples:
                    print(f"  {sample}: custom sample")
                    continue
                print(
                    f"  [WARNING] {sample} is not present in process_names.yaml "
                    "or process_groups.yaml"
                )
                continue

            if sample == plot_groups_cfg.get("other_group", {}).get("key", "OTHER"):
                print("  OTHER: background group (unmapped backgrounds)")
                continue

            if sample in plot_groups_cfg.get("background_groups", {}):
                members = get_group_members(
                    plot_groups_cfg["background_groups"][sample]
                )
                print(f"  {sample}: background group ({', '.join(members)})")
                continue

            info = classify_plot_sample(sample, process_cfg, plot_groups_cfg)

            if info is None:
                print(f"  [SKIP] {sample}")
            else:
                print(f"  {sample}: {info['type']}")

    # =====================================================
    # Read input ROOT files
    # =====================================================

    loaded_overrides = set()

    def read_inputs(primary_input):
        """Every requested process of one input directory: {process: info}."""
        input_processes = {}
        all_found_variables = set()

        scan_roots = [(primary_input, False)]
        for override_path in sample_input_overrides.values():
            override_dir = os.path.dirname(override_path)
            if all(root != override_dir for root, _ in scan_roots):
                scan_roots.append((override_dir, True))

        for scan_root, override_only in scan_roots:
          for indir, subdirs, infiles in os.walk(scan_root):

            for inFile in sorted(infiles):
                if not inFile.endswith(".root"):
                    continue

                full_path = os.path.abspath(os.path.join(indir, inFile))

                process_name = normalize_sample_name(inFile)

                override_path = sample_input_overrides.get(process_name)
                if override_path is not None:
                    if full_path != override_path:
                        continue
                    loaded_overrides.add(process_name)
                    print(f"[INFO] Alternate input for {process_name}: {full_path}")
                elif override_only:
                    continue

                if requested_samples is not None and process_name not in requested_samples:
                    continue

                sample_info = classify_plot_sample(process_name, process_cfg, plot_groups_cfg)
                if sample_info is None and args.allow_custom_samples:
                    sample_info = make_custom_sample_info(
                        process_name,
                        index=len(input_processes),
                    )

                if sample_info is None:
                    print(
                        f"[WARNING] Sample {process_name} not configured "
                        f"or skipped. File: {full_path}"
                    )
                    continue

                # print(
                #     f"[INFO] Loading {process_name} "
                #     f"as {sample_info['type']} from {full_path}"
                # )

                root_file = ROOT.TFile.Open(full_path, "READ")

                if not root_file or root_file.IsZombie():
                    print(f"[WARNING] Invalid ROOT file: {full_path}")
                    continue

                input_processes[process_name] = {
                    "input": full_path,
                    # The primary file may itself be a merged nominal+systematics
                    # file.  Keeping it as a systematic source also lets combined
                    # eras recover each decorrelated variation from the matching
                    # physical-era file below the same base directory.
                    "systematic_inputs": [full_path],
                    "color": sample_info["color"],
                    "name": sample_info["name"],
                    "is_data": sample_info["is_data"],
                    "is_signal": sample_info["is_signal"],
                    "type": sample_info["type"],
                    "hists": {
                        region_path: {},
                    },
                }

                scale_factor = process_scale_factors.get(process_name, 1.0)
                if args.dy_012j_weights != "none":
                    if root_file_has_patch(root_file):
                        print(
                            f"[WARNING] {full_path} already contains a persistent "
                            "DY component patch; the plot-time 0/1/2J weight is "
                            "not applied again."
                        )
                    else:
                        scale_factor *= dy_012j_plot_scale(
                            process_name, dy_012j_weights
                        )
                # print(process_name, scale_factor)

                available_hists = get_available_histograms(
                    root_file,
                    region_path,
                    scale_factor=scale_factor,
                    recursive=True,
                    exclude_2d=args.exclude_2d,
                ) 

                for available_hist, hist_name in available_hists:

                    base_name = hist_name
                    # JER component names intentionally do not carry the CMS_
                    # prefix (for example DNN_NNOutput_JEReta0pt02022Up).
                    # Treat them like every other shifted template so they inherit
                    # the nominal variable configuration instead of being skipped
                    # as unknown standalone variables.
                    for systematic_marker in ("_CMS_", "_QCD_", "_pdf_", "_JEReta"):
                        if systematic_marker in hist_name:
                            base_name = hist_name.split(systematic_marker, 1)[0]
                            break

                    if requested_variables_set is not None:
                        if hist_name not in requested_variables_set and base_name not in requested_variables_set:
                            continue

                    # Shifted histograms inherit binning and display settings from
                    # their nominal variable.  They are not separate entries in
                    # the histogram configuration.
                    try:
                        var_entry = findBinEntry(hist_cfg, base_name)
                    except KeyError:
                        print(
                            f"[WARNING] No configuration found "
                            f"for {hist_name}. Skip."
                        )
                        continue
                

                    if "x_rebin" in hist_cfg[var_entry]:
                        bins_to_compute = findNewBins(
                            hist_cfg,
                            var_entry,
                            dir_name=region_path,
                        )
                        new_bins = getNewBins(bins_to_compute)
                    else:
                        new_bins = hist_cfg[var_entry].get("x_bins", [])

                    rebinned_hist = available_hist

                    if args.rebin:
                        rebinned_hist = RebinHisto(
                            available_hist,
                            new_bins,
                            process_name,
                            # Keep under/overflow in their ROOT bins so they are
                            # available to the yield, but do not fold them into
                            # the first/last visible plotting bin.
                            wantOverflow=False,
                        )

                    if rebinned_hist is None:
                        continue

                    rebinned_hist.SetDirectory(0)

                    if is_valid_histogram(rebinned_hist):
                        input_processes[process_name]["hists"][region_path][hist_name] = (
                            rebinned_hist
                        )
                        # Systematic Up/Down shapes are stored so that the
                        # uncertainty machinery can retrieve them, but only the
                        # nominal variable is a standalone plot candidate.
                        all_found_variables.add(base_name)

                root_file.Close()
        return input_processes, all_found_variables

    if args.era_inputs:
        # Several eras summed here, without a merged file: see sum_era_inputs.
        era_dirs = [os.path.abspath(d) for d in args.era_inputs.split(",") if d]
        input_processes, all_found_variables = sum_era_inputs(
            [read_inputs(d) for d in era_dirs], region_path)
    else:
        input_processes, all_found_variables = read_inputs(os.path.abspath(args.input))

    for sample, path in sample_input_overrides.items():
        if sample not in loaded_overrides:
            print(f"[WARNING] Alternate input for {sample} was not loaded: {path}")

    # Load shifted templates from independent systematic directories.  This
    # avoids building one very large all-systematics ROOT file per process.
    for source_spec in args.systematic_input:
        if "=" not in source_spec:
            raise ValueError(
                f"--systematic-input expects NAME=DIR, got {source_spec!r}"
            )
        source_name, source_dir = source_spec.split("=", 1)
        print(f"[SYST INPUT] {source_name}: {source_dir}")
        for process_name, process_info in input_processes.items():
            source_path = os.path.join(source_dir, f"{process_name}.root")
            if not os.path.isfile(source_path):
                continue
            process_info["systematic_inputs"].append(source_path)
            source_file = ROOT.TFile.Open(source_path, "READ")
            if not source_file or source_file.IsZombie():
                print(f"[WARNING] Invalid systematic ROOT file: {source_path}")
                continue
            available_hists = get_available_histograms(
                source_file,
                region_path,
                scale_factor=process_scale_factors.get(process_name, 1.0),
                recursive=True,
                exclude_2d=args.exclude_2d,
            )
            for available_hist, hist_name in available_hists:
                if not hist_name.endswith(("Up", "Down")):
                    continue
                base_name = re.sub(
                    r"_(?:CMS_|QCD_|pdf_|JEReta).*(?:Up|Down)$", "", hist_name
                )
                if requested_variables_set is not None and base_name not in requested_variables_set:
                    continue
                try:
                    var_entry = findBinEntry(hist_cfg, base_name)
                except KeyError:
                    continue
                rebinned_hist = available_hist
                if args.rebin:
                    if "x_rebin" in hist_cfg[var_entry]:
                        bins_to_compute = findNewBins(
                            hist_cfg, var_entry, dir_name=region_path
                        )
                        new_bins = getNewBins(bins_to_compute)
                    else:
                        new_bins = hist_cfg[var_entry].get("x_bins", [])
                    rebinned_hist = RebinHisto(
                        available_hist,
                        new_bins,
                        process_name,
                        wantOverflow=False,
                    )
                if rebinned_hist is None:
                    continue
                rebinned_hist.SetDirectory(0)
                if is_valid_histogram(rebinned_hist):
                    process_info["hists"][region_path][hist_name] = rebinned_hist
            source_file.Close()

    input_processes = apply_plot_groups(
        input_processes,
        plot_groups_cfg,
        active_group_names=requested_plot_groups,
    )

    for process_name, process_info in input_processes.items():
        if process_name in sample_color_overrides:
            process_info["color"] = sample_color_overrides[process_name]
        if process_name in sample_label_overrides:
            process_info["name"] = sample_label_overrides[process_name]

    # =====================================================
    # Summary
    # =====================================================

    # print("\nLoaded samples:")

    # for sample, info in input_processes.items():
    #     n_hists = len(info["hists"].get(region_path, {}))
    #     print(
    #         f"  {sample:20s} "
    #         f"type={info['type']:10s} "
    #         f"n_hists={n_hists}"
    #     )

    # =====================================================
    # Produce plots
    # =====================================================

    if len(all_found_variables) == 0:

        if requested_variables is not None:
            print(
                "[ERROR] None of the requested variables were found: "
                f"{', '.join(requested_variables)}"
            )

        print(
            f"[ERROR] No valid histogram found "
            f"for region {region_path}."
        )

    else:

        output_dir_path = os.path.join(
            args.output,
            args.era,
            region_path,
        )

        os.makedirs(output_dir_path, exist_ok=True)

        variables_to_plot = sorted(all_found_variables)

        if requested_variables is not None:
            missing_variables = [
                variable
                for variable in requested_variables
                if variable not in all_found_variables
            ]

            if missing_variables:
                print(
                    "[WARNING] Requested variables not found in "
                    f"{region_path}: {', '.join(missing_variables)}"
                )

            variables_to_plot = [
                variable
                for variable in requested_variables
                if variable in all_found_variables
            ]

        print(
            f"\n--> Generazione di {len(variables_to_plot)} "
            f"plot strutturati in corso..."
        )

        multipage_pdf = None
        if args.multipage_pdf is not None:
            multipage_path = args.multipage_pdf
            if not os.path.isabs(multipage_path):
                multipage_path = os.path.join(output_dir_path, multipage_path)
            if not multipage_path.lower().endswith(".pdf"):
                multipage_path += ".pdf"
            os.makedirs(os.path.dirname(multipage_path), exist_ok=True)
            multipage_pdf = PdfPages(multipage_path)
            print(f"[multipage PDF] {multipage_path}")

        try:
            for variable in variables_to_plot:

                plot_base_path = os.path.join(
                    output_dir_path,
                    variable,
                )

                os.makedirs(
                    os.path.dirname(plot_base_path),
                    exist_ok=True,
                )

                make_stacked_plot(
                    samples_dict=input_processes,
                    config_page=config_setup,
                    category=region_path,
                    variable=variable,
                    out_name=plot_base_path,
                    want_data=args.wantData,
                    do_stack=args.do_stack,
                    fill_hists=args.fill_hists,
                    ratio_reference=args.ratio_reference,
                    normalize_dy_to_data=args.normalize_dy_to_data,
                    normalize_mc_to_data=args.normalize_mc_to_data,
                    era=args.era,
                    dy_normalization_sample=args.dy_normalization_sample,
                    dy_composition=args.dy_composition,
                    dy_component_reweighted=(args.dy_012j_weights != "none"),
                    show_systematics=args.systematics,
                    systematic_groups=args.systematicGroup,
                    combined_eras=combined_eras,
                    overlay_systematic=args.overlaySystematic,
                    log_uncertainties=args.logUncertainties,
                    include_total_systematics=args.totalSystematics,
                    show_mc_stat_uncertainty=not args.noMCStatUncertainty,
                    multipage_pdf=multipage_pdf,
                )
        finally:
            if multipage_pdf is not None:
                multipage_pdf.close()

    print(
        f"\n[SUCCESS] Processing completed "
        f"in {time.time() - startTime:.2f} seconds."
    )
