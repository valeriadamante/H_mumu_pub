#!/usr/bin/env python3
import argparse
import copy
import gc
import os
import shutil
import sys
import time
import traceback
from pathlib import Path
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.EnableThreadSafety()
sys.path.append(os.environ["ANALYSIS_PATH"])
import common.utilities as utilities
from common.add_vars import GetSelectionSuffixForSystematic
from common.prepare_rdf import prepare_rdf, prepare_region_dataframes
from common.dnn_histogram_production import (
    needs_sideband_mass_shift,
    shifted_output_column,
)
from common.jet_component_splitting import (
    DY_COMPONENT_FILE_LABELS,
    GGF_COMPONENT_VARIABLES,
    VBF_COMPONENTS,
    VBF_ETA_REGIONS,
    add_jet_component_categories,
    add_vbf_eta_region_categories,
    component_output_directory,
    expanded_jet_component_categories,
    jet_components_enabled_for_dataset,
    variable_for_component,
)
from common.utilities import (
    read_manifest,
    initialize_root_runtime,
    validate_file,
    GetModel,
    findBinEntry,
    list_root_files,
    get_segmentation_dict,
)
from common.systematic_correlations import nuisance_name as correlated_nuisance_name
from corrections.qcd_scale import get_qcd_scale_points
initialize_root_runtime()
_METADATA_CACHE = {}
def profile_log(job, phase, started_at):
    elapsed = time.perf_counter() - started_at
    print(f"[PROFILE][{job}] {phase}: {elapsed:.3f} s", flush=True)
    return time.perf_counter()
def parse_file_shard(value):
    """Parse the 'INDEX/TOTAL' form of --file-shard into a 1-based pair."""
    parts = str(value).split("/")
    if len(parts) != 2:
        raise ValueError("--file-shard must be given as INDEX/TOTAL")
    try:
        shard_index, shard_total = (int(part) for part in parts)
    except ValueError:
        raise ValueError("--file-shard must be given as INDEX/TOTAL") from None
    if shard_total < 1:
        raise ValueError("--file-shard TOTAL must be >= 1")
    if not 1 <= shard_index <= shard_total:
        raise ValueError("--file-shard INDEX must satisfy 1 <= INDEX <= TOTAL")
    return shard_index, shard_total
def batch_dict(items, batch_size):
    entries = list(items.items())
    return [
        dict(entries[index : index + batch_size])
        for index in range(0, len(entries), batch_size)
    ]
def batch_list(items, batch_size):
    return [
        items[index : index + batch_size]
        for index in range(0, len(items), batch_size)
    ]
def safe_mkdir(path):
    if path:
        os.makedirs(path, exist_ok=True)
def remove_file_if_exists(path):
    if path and os.path.exists(path):
        try:
            os.remove(path)
        except Exception as e:
            print(f"[WARNING] Could not remove file {path}: {e}")
def histogram_directory_path(mass_region, category):
    if category.startswith("VBF_eta_"):
        eta_region = category.removeprefix("VBF_eta_")
        return f"{mass_region}_VBF/{eta_region}"
    return f"{mass_region}_{category}"
def copy_root_directory(source_file, source_path, target_file, target_path):
    """Copy all objects from one ROOT directory into a nested target directory."""
    source = source_file.GetDirectory(source_path)
    if not source:
        return False
    target = utilities.mkdir_recursive(target_file, target_path)
    for key in source.GetListOfKeys():
        obj = key.ReadObj()
        target.cd()
        target.WriteTObject(obj, obj.GetName(), "Overwrite")
    return True
def assert_complete_root_file(path):
    """Raise if path is not a complete ROOT file (zombie, recovered or without keys)."""
    handle = ROOT.TFile.Open(str(path), "READ")
    try:
        if (not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered)
                or handle.GetListOfKeys().GetSize() == 0):
            raise RuntimeError(f"Incomplete ROOT output (zombie, recovered or empty): {path}")
    finally:
        if handle:
            handle.Close()


def split_dy_jet_component_outputs(
    output_file,
    mass_regions,
    process_label="DY",
    include_vbf_eta_regions=False,
    requested_categories=None,
):
    """Create inclusive and per-component PU/hard files with fit-ready layout."""
    requested = set(requested_categories or ("ggF", "VBF"))
    passthrough_categories = requested - {"ggF", "VBF"}
    eta_regions = VBF_ETA_REGIONS if include_vbf_eta_regions else ("incl",)
    output_path = Path(output_file)
    source = ROOT.TFile.Open(str(output_path), "READ")
    # A recovered source has lost keys: the split would rewrite it as a clean but incomplete file.
    if not source or source.IsZombie() or source.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"Could not open DY component staging file (zombie or recovered): {output_file}")
    inclusive_tmp = output_path.with_name(f".{output_path.name}.inclusive.tmp.root")
    inclusive = ROOT.TFile.Open(str(inclusive_tmp), "RECREATE")
    component_files = {}
    files_by_label = {}
    try:
        selected_components = [
            component
            for component in DY_COMPONENT_FILE_LABELS
            if (component.startswith("ggF_") and "ggF" in requested)
            or (component.startswith("VBF_") and "VBF" in requested)
        ]
        for component in selected_components:
            label = DY_COMPONENT_FILE_LABELS[component]
            if process_label != "DY" and label.startswith("DY_"):
                label = f"{process_label}_{label[len('DY_'):]}"
            if label not in files_by_label:
                component_path = output_path.with_name(
                    f"{output_path.stem}_{label}{output_path.suffix}"
                )
                files_by_label[label] = (
                    component_path,
                    ROOT.TFile.Open(str(component_path), "RECREATE"),
                )
            component_files[component] = files_by_label[label]
        for mass_region in mass_regions:
            for category in sorted(passthrough_categories):
                copy_root_directory(
                    source,
                    histogram_directory_path(mass_region, category),
                    inclusive,
                    histogram_directory_path(mass_region, category),
                )
            if "ggF" in requested:
                copy_root_directory(
                    source,
                    f"{mass_region}_DY_inclusive_ggF",
                    inclusive,
                    component_output_directory(mass_region, "ggF"),
                )
            if "VBF" in requested:
                for eta_region in eta_regions:
                    copy_root_directory(
                        source,
                        f"{mass_region}_DY_inclusive_VBF_{eta_region}",
                        inclusive,
                        component_output_directory(
                            mass_region,
                            "VBF",
                            eta_region if include_vbf_eta_regions else None,
                        ),
                    )
            if "ggF" in requested:
                for component in GGF_COMPONENT_VARIABLES:
                    _, target = component_files[component]
                    copy_root_directory(
                        source,
                        f"{mass_region}_{component}",
                        target,
                        component_output_directory(mass_region, "ggF"),
                    )
            if "VBF" in requested:
                for component in VBF_COMPONENTS:
                    _, target = component_files[component]
                    for eta_region in eta_regions:
                        copy_root_directory(
                            source,
                            f"{mass_region}_{component}_{eta_region}",
                            target,
                            component_output_directory(
                                mass_region,
                                "VBF",
                                eta_region if include_vbf_eta_regions else None,
                            ),
                        )
    finally:
        source.Close()
        inclusive.Close()
        for _, target in files_by_label.values():
            target.Close()
    os.replace(inclusive_tmp, output_path)
    print("[INFO] DY inclusive/component outputs:")
    print(f"[INFO]   inclusive: {output_path}")
    for component in selected_components:
        print(f"[INFO]   {component}: {component_files[component][0]}")
    return [path for path, _ in files_by_label.values()]
def format_systematic_info(syst_info, scale=None):
    formatted = {}
    for key, value in syst_info.items():
        if isinstance(value, str) and scale is not None:
            formatted[key] = value.replace("{scale}", scale)
        else:
            formatted[key] = value
    return formatted
def nuisance_histogram_name(variable, syst_name, syst_info, era, process):
    if syst_name == "Central":
        return variable
    nuisance_name = correlated_nuisance_name(
        dict(syst_info, name=syst_info.get("name", syst_name)), era,
        process=process,
        pdf_process=pdf_process_label(syst_info.get("pdf_config", {}), process),
    )
    direction = syst_info.get("direction")
    if direction:
        nuisance_name = f"{nuisance_name}{direction.capitalize()}"
    return f"{variable}_{nuisance_name}"
def unique_metadata_inputs(paths):
    unique_paths = []
    seen = set()
    for path in paths:
        if not path:
            continue
        normalized = os.path.abspath(path)
        if normalized in seen:
            continue
        seen.add(normalized)
        unique_paths.append(path)
    return unique_paths
def expand_metadata_inputs(paths):
    expanded = []
    for raw_path in unique_metadata_inputs(paths):
        path = Path(raw_path)
        if path.is_dir():
            expanded.extend(str(candidate) for candidate in sorted(path.rglob("*.json")))
        elif path.is_file():
            expanded.append(str(path))
        else:
            raise FileNotFoundError(f"Metadata input does not exist: {raw_path}")
    return unique_metadata_inputs(expanded)
def get_combined_segmentation_dict(
    input_paths,
    node="gen",
    fallback_to_initial=True,
    warn_if_missing=True,
):
    metadata_inputs = unique_metadata_inputs(input_paths)
    cache_key = (
        tuple(os.path.abspath(path) for path in metadata_inputs),
        node,
        bool(fallback_to_initial),
    )
    if cache_key in _METADATA_CACHE:
        return _METADATA_CACHE[cache_key]
    combined = get_segmentation_dict(
        metadata_inputs,
        node=node,
        fallback_to_initial=fallback_to_initial,
        warn_if_missing=warn_if_missing,
    )
    _METADATA_CACHE[cache_key] = combined
    return combined
def get_qcd_scale_segmentation_dict(input_paths, point_name, **kwargs):
    warn_if_missing = kwargs.pop("warn_if_missing", True)
    node_candidates = [
        f"qcd_scale__{point_name}",
        f"gen_qcdScale_{point_name}",
    ]
    for node in node_candidates:
        sums = get_combined_segmentation_dict(
            input_paths,
            node=node,
            warn_if_missing=False,
            **kwargs,
        )
        if sums:
            return sums
    if warn_if_missing:
        print(
            "[WARNING] No QCD scale segmentation JSON information found for "
            f"{point_name} under: " + ", ".join(unique_metadata_inputs(input_paths))
        )
    return {}
def get_qcd_scale_variations(qcd_scale_config):
    return qcd_scale_config.get(
        "variations",
        [
            {
                "name": "QCDscaleMuR_{process}",
                "down": "muR0p5_muF1",
                "up": "muR2_muF1",
            },
            {
                "name": "QCDscaleMuF_{process}",
                "down": "muR1_muF0p5",
                "up": "muR1_muF2",
            },
        ],
    )
def get_qcd_scale_source_points(qcd_scale_config):
    points_by_name = {
        point["name"]: point
        for point in get_qcd_scale_points(qcd_scale_config)
    }
    source_names = []
    for variation in get_qcd_scale_variations(qcd_scale_config):
        for direction in ("down", "up"):
            point_name = variation[direction]
            if point_name not in points_by_name:
                raise ValueError(
                    f"QCD scale variation point '{point_name}' is not defined "
                    "in qcd_scale.points."
                )
            if point_name not in source_names:
                source_names.append(point_name)
    return [points_by_name[name] for name in source_names]
def qcd_scale_process_label(qcd_scale_config, process):
    process_labels = qcd_scale_config.get("process_labels", {})
    if process in process_labels:
        return process_labels[process]
    lower_process = process.lower()
    if lower_process.startswith(("dy", "w")) or "ewk" in lower_process:
        return "V"
    if lower_process.startswith(("tt", "st", "tw")):
        return "ttbar"
    if lower_process.startswith("vbf") or "vbfh" in lower_process:
        return "qqH"
    if lower_process.startswith(("gluglu", "ggh")):
        return "ggH"
    if lower_process.startswith("vh") or "zh" in lower_process or "wh" in lower_process:
        return "VH"
    if lower_process.startswith(("tth", "ttH".lower())):
        return "ttH"
    if lower_process.startswith("vvv"):
        return "VVV"
    if lower_process.startswith("vv") or lower_process in {"ww", "wz", "zz"}:
        return "VV"
    return process
def pdf_process_label(pdf_config, process):
    process_labels = pdf_config.get("process_labels", {})
    if process in process_labels:
        return process_labels[process]
    lower_process = process.lower()
    if lower_process.startswith(("dy", "w")) or "ewk" in lower_process:
        return "qqbar"
    if lower_process.startswith("vbf") or "vbfh" in lower_process:
        return "Higgs_qqH"
    if lower_process.startswith(("gluglu", "ggh")):
        return "Higgs_ggH"
    if lower_process.startswith("vh") or "zh" in lower_process or "wh" in lower_process:
        return "Higgs_VH"
    if lower_process.startswith("tth"):
        return "Higgs_ttH"
    if lower_process.startswith("tt"):
        return "gg"
    if lower_process.startswith(("st", "tw")):
        return "gq"
    if lower_process.startswith("vvv"):
        return "qqbar"
    if lower_process.startswith("vv") or lower_process in {"ww", "wz", "zz"}:
        return "qqbar"
    return process
def configure_available_qcd_scale(syst_cfg, input_paths, is_data, mode):
    qcd_config = syst_cfg.get("qcd_scale", {})
    if (
        is_data
        or mode == "central"
        or not qcd_config.get("enabled", False)
    ):
        return syst_cfg
    available_points = []
    missing_points = []
    source_points = get_qcd_scale_source_points(qcd_config)
    for point in source_points:
        point_name = point["name"]
        sums = get_qcd_scale_segmentation_dict(
            input_paths,
            point_name,
            fallback_to_initial=False,
            warn_if_missing=False,
        )
        if sums:
            available_points.append(point)
        else:
            missing_points.append(point_name)
    configured = copy.deepcopy(syst_cfg)
    configured["qcd_scale"]["points"] = available_points
    available_names = {point["name"] for point in available_points}
    configured["qcd_scale"]["variations"] = [
        variation
        for variation in get_qcd_scale_variations(qcd_config)
        if variation["down"] in available_names and variation["up"] in available_names
    ]
    if not missing_points:
        return configured
    missing_policy = qcd_config.get("missing_sums", "error")
    message = (
        "QCD scale sums are missing from the skim reports for: "
        + ", ".join(missing_points)
    )
    if missing_policy == "error":
        raise RuntimeError(
            f"{message}. Reproduce the skim with --want-variations."
        )
    if missing_policy != "skip":
        raise ValueError(
            "qcd_scale.missing_sums must be either 'skip' or 'error'"
        )
    if available_points:
        print(
            f"[WARNING] {message}. Missing QCD scale points will be skipped; "
            f"{len(configured['qcd_scale']['variations'])} QCD scale "
            "variation(s) will still be produced."
        )
    else:
        configured["qcd_scale"]["enabled"] = False
        print(
            f"[WARNING] {message}. QCD scale templates will be skipped; "
            "all other requested systematics will still be produced."
        )
    return configured
def get_systs_to_run(syst_cfg, mode):
    systs_to_run = {
        "Central": syst_cfg["systematics"]["Central"]
    }
    if mode == "central":
        return systs_to_run
    scales = syst_cfg.get("scales", ["up", "down"])
    for syst_name, syst_info in syst_cfg.get("systematics", {}).items():
        if syst_name == "Central":
            continue
        if mode == "jec-jer" and syst_name not in ("JER", "JES_Total"):
            continue
        components = syst_info.get("components", ())
        for component in components or (None,):
            for scale in scales:
                output_base = component or syst_name
                output_name = f"{output_base}{scale.capitalize()}"
                formatted = format_systematic_info(syst_info, scale=scale)
                formatted.pop("components", None)
                if component:
                    formatted["jer_component"] = component
                    formatted["source_jet_suffix"] = formatted["jet_suffix"]
                    formatted["jet_suffix"] = f"_{component}{scale.capitalize()}"
                    formatted["name"] = f"{component}{{era}}"
                formatted["direction"] = scale
                systs_to_run[output_name] = formatted
    for weight_name, weight_info in syst_cfg.get("weights", {}).items():
        if weight_name == "Central":
            continue
        if weight_info.get("derived_envelope", False):
            continue
        if "{scale}" in weight_name:
            for scale in scales:
                output_name = weight_name.format(scale=scale)
                formatted = format_systematic_info(weight_info, scale=scale)
                formatted["direction"] = scale
                if weight_name.startswith("PDF_"):
                    formatted["pdf_config"] = syst_cfg.get("pdf", {})
                systs_to_run[output_name] = formatted
        else:
            formatted = dict(weight_info)
            for direction in scales:
                if weight_name.endswith(f"_{direction}"):
                    formatted["direction"] = direction
                    break
            if weight_name.startswith("PDF_"):
                formatted["pdf_config"] = syst_cfg.get("pdf", {})
            systs_to_run[weight_name] = formatted
    qcd_scale_config = syst_cfg.get("qcd_scale", {})
    if qcd_scale_config.get("enabled", False):
        for point in get_qcd_scale_source_points(qcd_scale_config):
            point_name = point["name"]
            output_name = f"QCDScale__{point_name}"
            systs_to_run[output_name] = {
                "jet_suffix": "",
                "muon_suffix": "",
                "name": output_name,
                "weight": f"weight__{output_name}",
            }
    return systs_to_run
def parse_requested_systematics(values):
    if not values:
        return []
    requested = []
    for value in values:
        for item in str(value).replace(",", " ").split():
            if item:
                requested.append(item)
    return list(dict.fromkeys(requested))
def expand_systematic_group_alias(requested_name, available_systematics):
    normalized_name = requested_name.lower().replace("_", "").replace("-", "")
    aliases = {
        "jesregrouped": tuple(name for name in available_systematics
                              if name.startswith(("JESRegrouped_", "JESRelativeSample_"))),
        "jerc": (
            *tuple(name for name in available_systematics if name.startswith("JER")),
            "JES_TotalUp", "JES_TotalDown",
        ),
        "qcdscale": tuple(
            name
            for name in available_systematics
            if name.startswith("QCDScale__")
        ),
        "qcdscales": tuple(
            name
            for name in available_systematics
            if name.startswith("QCDScale__")
        ),
        "pdf": tuple(
            name
            for name in available_systematics
            if name.startswith("PDF_")
        ),
        "scaleweight": tuple(
            name
            for name in available_systematics
            if name.startswith("QCDScale__") or name.startswith("PDF_")
        ),
        "scare": (
            "MuonScaleUp",
            "MuonScaleDown",
            "MuonResUp",
            "MuonResDown",
        ),
        "muon": (
            "MuonID_up",
            "MuonID_down",
            "MuonIso_up",
            "MuonIso_down",
            "singleMuTrigger_up",
            "singleMuTrigger_down",
        ),
        "pu": (
            "PU_up",
            "PU_down",
        ),
    }
    expanded = [
        name
        for name in aliases.get(normalized_name, ())
        if name in available_systematics
    ]
    if expanded or requested_name in available_systematics:
        return expanded

    # Campaigns are grouped by nuisance family, while the histogram config
    # stores the concrete Up/Down variations.
    return [
        name
        for name in (f"{requested_name}Up", f"{requested_name}Down")
        if name in available_systematics
    ]
def filter_systs_to_run(systs_to_run, requested_systematics):
    requested_raw = parse_requested_systematics(requested_systematics)
    requested = []
    for name in requested_raw:
        expanded = expand_systematic_group_alias(name, systs_to_run)
        if expanded:
            requested.extend(expanded)
        else:
            requested.append(name)
    requested = list(dict.fromkeys(requested))
    if not requested:
        return systs_to_run
    missing = [name for name in requested if name not in systs_to_run]
    if missing:
        available = ", ".join(sorted(systs_to_run))
        raise ValueError(
            "Unknown requested systematic(s): "
            + ", ".join(missing)
            + ". Available systematics: "
            + available
        )
    return {name: systs_to_run[name] for name in requested}
def validate_systematic_isolation(systs_to_run):
    """Prevent one nuisance family from shifting unrelated inputs."""
    for name, info in systs_to_run.items():
        jet_suffix = info.get("jet_suffix", "")
        muon_suffix = info.get("muon_suffix", "")
        weight = info.get("weight", "weight__Central")
        if name.startswith(("JER", "JES")) and (
            muon_suffix or weight != "weight__Central"
        ):
            raise ValueError(
                f"{name} must vary only jet_suffix (found muon_suffix="
                f"{muon_suffix!r}, weight={weight!r})"
            )
        if name.startswith(("MuonScale", "MuonRes")) and (
            jet_suffix or weight != "weight__Central"
        ):
            raise ValueError(
                f"{name} must vary only muon_suffix (found jet_suffix="
                f"{jet_suffix!r}, weight={weight!r})"
            )
        if name.startswith(("MuonID_", "MuonIso_", "singleMuTrigger_")) and (
            jet_suffix or muon_suffix
        ):
            raise ValueError(
                f"{name} must vary only its weight (found jet_suffix="
                f"{jet_suffix!r}, muon_suffix={muon_suffix!r})"
            )
def write_qcd_scale_variations(
    output_file,
    syst_cfg,
    variables,
    mass_regions,
    categories,
    era,
    process,
):
    qcd_scale_config = syst_cfg.get("qcd_scale", {})
    if not qcd_scale_config.get("enabled", False):
        return
    variations = get_qcd_scale_variations(qcd_scale_config)
    process_label = qcd_scale_process_label(qcd_scale_config, process)
    source_suffixes = sorted(
        {
            f"QCDScale__{variation[direction]}"
            for variation in variations
            for direction in ("down", "up")
        }
    )
    for mass_region in mass_regions:
        for category in categories:
            directory = output_file.GetDirectory(
                histogram_directory_path(mass_region, category)
            )
            if not directory:
                continue
            for variable in variables:
                for variation in variations:
                    nuisance_name = correlated_nuisance_name(
                        dict(qcd_scale_config, **variation), era,
                        process=process_label,
                    )
                    for direction, shape_direction in (
                        ("down", "Down"),
                        ("up", "Up"),
                    ):
                        source = directory.Get(
                            f"{variable}_QCDScale__{variation[direction]}"
                        )
                        if not source or not source.InheritsFrom(ROOT.TH1.Class()):
                            continue
                        hist = source.Clone(
                            f"{variable}_{nuisance_name}{shape_direction}"
                        )
                        hist.SetDirectory(0)
                        directory.WriteTObject(hist, hist.GetName(), "Overwrite")
                for source_suffix in source_suffixes:
                    directory.Delete(
                        f"{variable}_{source_suffix};*"
                    )
def get_histogram_variable(variable, syst_info, available_columns):
    jet_suffix = syst_info.get("jet_suffix", "")
    muon_suffix = syst_info.get("muon_suffix", "")
    candidates = []
    if jet_suffix:
        candidates.append(f"{variable}{jet_suffix}")
    if muon_suffix:
        candidates.append(f"{variable}{muon_suffix}")
    candidates.append(variable)
    return next(
        (candidate for candidate in candidates if candidate in available_columns),
        None,
    )
def produce_histograms(args_tuple):
    (
        input_files,
        dataset_seg_dict,
        dataset_qcd_scale_seg_dicts,
        args,
        is_data,
        sel_cfg,
        syst_cfg,
        vars_to_make_hist,
        masses_regions,
        masses_regions_list,
        categories,
        categories_list,
        hist_cfg,
        systs_to_run,
        dnn_payloads,
        btag_algo,
    ) = args_tuple
    output_path = args.output_file
    out_file = None
    try:
        if args.rdf_threads > 1 and not ROOT.IsImplicitMTEnabled():
            ROOT.EnableImplicitMT(args.rdf_threads)
        print(
            f"[JOB {args.dataset_name}] Starting with {len(input_files)} file(s)"
        )
        phase_started = time.perf_counter()
        profile_log(args.dataset_name, "input setup", phase_started)
        seg_dict = None if is_data else dataset_seg_dict
        qcd_scale_seg_dicts = (
            {} if is_data else dataset_qcd_scale_seg_dicts
        )
        if is_data:
            print(
                f"[JOB {args.dataset_name}] Data dataset: "
                "segmentation metadata disabled"
            )
        else:
            print(
                f"[JOB {args.dataset_name}] Using "
                f"{len(seg_dict)} segmentation entries for "
                f"{len(input_files)} ROOT file(s)"
            )
        rdf_started = time.perf_counter()
        prepared = prepare_rdf(
            dataset_name=args.dataset_name, era=args.era,
            selections_cfg=sel_cfg, systematics_cfg=syst_cfg,
            is_data=is_data, input_dir=args.root_input, input_files=input_files,
            seg_dict=seg_dict, qcd_scale_seg_dicts=qcd_scale_seg_dicts,
            systs_to_run=systs_to_run,
            want_variations=args.systematics_mode != "central",
            dnn_payloads=dnn_payloads, btag_algo=btag_algo,
            additional_cuts=args.additional_cuts, dnn_model_set=args.dnn_model_set,
            skip_validation=True,
            split_jet_multiplicity=args.dy_jet_components,
            include_vbf_eta_regions=args.vbf_eta_regions,
            component_categories=getattr(args, "pu_hard_requested_categories", ("ggF", "VBF")),
            vbf_pair_definition=args.vbf_pair_definition,
        )
        rdf_base = prepared.get("inclusive")
        profile_log(args.dataset_name, "RDataFrame preparation", rdf_started)
        booking_setup_started = time.perf_counter()
        from common.utilities import dataset_region_allowed
        stored_regions = [
            name
            for name, info in masses_regions.items()
            if name in masses_regions_list and info.get("store", False)
            and (not args.region_sample_routing or dataset_region_allowed(args.dataset_name, name))
        ]
        stored_categories = [
            name
            for name, info in categories.items()
            if name in categories_list and info.get("store", False)
        ]
        hist_specs = {}
        for variable in vars_to_make_hist:
            config_key = findBinEntry(hist_cfg, variable)
            configured_columns = hist_cfg[config_key].get("var_list")
            columns = tuple(configured_columns or (variable,))
            hist_specs[variable] = {
                "columns": columns,
                "model": GetModel(
                    hist_cfg, variable, dims=len(columns), era=args.era
                ),
            }
        base_columns = (
            {str(column) for column in rdf_base.GetColumnNames()}
            if rdf_base is not None
            else set()
        )
        filtered_rdfs = prepare_region_dataframes(
            rdf_base, mass_regions=stored_regions, categories=stored_categories,
            systs_to_run=systs_to_run, variables=vars_to_make_hist,
            btag_algo=btag_algo, era=args.era, dnn_model_set=args.dnn_model_set,
        )
        profile_log(args.dataset_name, "selection and DNN graph construction", booking_setup_started)
        output_open_started = time.perf_counter()
        out_file = ROOT.TFile(output_path, "RECREATE")
        if not out_file or out_file.IsZombie():
            raise RuntimeError(f"Could not create output file: {output_path}")
        directories = {
            (mass_region, category): utilities.mkdir_recursive(
                out_file, histogram_directory_path(mass_region, category)
            )
            for mass_region in stored_regions
            for category in stored_categories
        }
        profile_log(args.dataset_name, "output open/directory creation", output_open_started)
        systematic_batches = batch_dict(systs_to_run, args.systematic_batch_size)
        variable_batches = batch_list(
            list(hist_specs), args.variable_batch_size
        )
        work_batches = [
            (systematic_batch, variable_batch)
            for systematic_batch in systematic_batches
            for variable_batch in variable_batches
        ]
        total_booked = 0
        for batch_index, (systematic_batch, variable_batch) in enumerate(
            work_batches, start=1
        ):
            histogram_booking_started = time.perf_counter()
            booked_hists = []
            for syst_name, syst_info in systematic_batch.items():
                weight_name = syst_info["weight"]
                selection_suffix = GetSelectionSuffixForSystematic(
                    syst_name, syst_info
                )
                if rdf_base is not None and weight_name not in base_columns:
                    raise RuntimeError(
                        f"Weight column '{weight_name}' not found for systematic "
                        f"'{syst_name}'"
                    )
                for mass_region in stored_regions:
                    for category in stored_categories:
                        rdf_filtered = None
                        available_columns = base_columns
                        filtered_entry = filtered_rdfs.get(
                            (mass_region, category, selection_suffix)
                        )
                        if filtered_entry is not None:
                            rdf_filtered, available_columns = filtered_entry
                        directory = directories[(mass_region, category)]
                        category_variables = set(
                            variable_for_component(
                                category, args.vbf_component_variables
                            )
                            if args.dy_jet_components
                            else vars_to_make_hist
                        )
                        for variable in variable_batch:
                            spec = hist_specs[variable]
                            if variable not in category_variables:
                                continue
                            model = spec["model"]
                            hist_name = nuisance_histogram_name(
                                variable,
                                syst_name,
                                syst_info,
                                args.era,
                                args.process_name,
                            )
                            hist_columns = tuple(
                                get_histogram_variable(
                                    column, syst_info, available_columns
                                )
                                for column in spec["columns"]
                            )
                            if needs_sideband_mass_shift(mass_region, variable):
                                hist_columns = (shifted_output_column(mass_region),)
                            if (
                                rdf_filtered is not None
                                and all(
                                    column in available_columns
                                    for column in hist_columns
                                )
                            ):
                                if len(hist_columns) == 1:
                                    hist_ptr = rdf_filtered.Histo1D(
                                        model, hist_columns[0], weight_name
                                    )
                                elif len(hist_columns) == 2:
                                    hist_ptr = rdf_filtered.Histo2D(
                                        model,
                                        hist_columns[0],
                                        hist_columns[1],
                                        weight_name,
                                    )
                                else:
                                    raise RuntimeError(
                                        f"Unsupported histogram dimension for "
                                        f"{variable}: {len(hist_columns)}"
                                    )
                                booked_hists.append(
                                    (directory, hist_name, hist_ptr, True)
                                )
                                continue
                            if rdf_filtered is not None:
                                print(
                                    f"[JOB {args.dataset_name}] WARNING: variable "
                                    f"'{variable}' not found for systematic "
                                    f"'{syst_name}'. Booking empty histogram."
                                )
                            hist = model.GetHistogram().Clone(hist_name)
                            hist.SetTitle(hist_name)
                            hist.Reset("ICES")
                            hist.SetDirectory(0)
                            booked_hists.append(
                                (directory, hist_name, hist, False)
                            )
            total_booked += len(booked_hists)
            batch_label = f"{batch_index}/{len(work_batches)}"
            profile_log(
                args.dataset_name,
                f"histogram booking batch {batch_label} "
                f"({len(booked_hists)} histograms)",
                histogram_booking_started,
            )
            print(
                f"[JOB {args.dataset_name}] Batch {batch_label}: "
                f"booked {len(booked_hists)} histograms for "
                f"{len(systematic_batch)} systematic variations and "
                f"{len(variable_batch)} variables."
            )
            if not booked_hists:
                continue
            histogram_actions = [
                hist_obj
                for _, _, hist_obj, needs_getvalue in booked_hists
                if needs_getvalue
            ]
            event_loop_started = time.perf_counter()
            if histogram_actions:
                ROOT.RDF.RunGraphs(histogram_actions)
            profile_log(
                args.dataset_name,
                f"ROOT event loop batch {batch_label} (RunGraphs)",
                event_loop_started,
            )
            output_write_started = time.perf_counter()
            for directory, hist_name, hist_obj, needs_getvalue in booked_hists:
                hist = hist_obj.GetValue() if needs_getvalue else hist_obj
                hist.SetName(hist_name)
                hist.SetTitle(hist_name)
                hist.SetDirectory(0)
                directory.cd()
                directory.WriteTObject(hist, hist_name, "Overwrite")
            out_file.Flush()
            profile_log(
                args.dataset_name,
                f"histogram materialization/write batch {batch_label}",
                output_write_started,
            )
            del histogram_actions
            del booked_hists
            gc.collect()
        if total_booked == 0:
            raise RuntimeError(
                "No histograms were booked. Check that --mass-regions and "
                "--categories are passed as space-separated values and match "
                "the selection configuration."
            )
        out_file.Close()
        out_file = None
        print(f"[JOB {args.dataset_name}] Done -> {output_path}")
        return output_path
    except Exception:
        if out_file:
            out_file.Close()
        traceback.print_exc()
        remove_file_if_exists(output_path)
        raise
    finally:
        # ApplyDNN stores predictions in a process-global C++ registry.  A
        # dataset job handles MC files serially, so retaining completed-file
        # payloads makes RSS grow monotonically until the cgroup kills it.
        from common.dnn_application import clear_prediction_registry
        clear_prediction_registry()
def main(argv=None):
    parser = argparse.ArgumentParser(description="Produce histograms from validated skimmed ROOT ntuples.")
    parser.add_argument("--era", required=True, help="Era, e.g. Run3_2022EE")
    parser.add_argument(
        "--root-input",
        "--input",
        dest="root_input",
        default=None,
        help="Skimmed ROOT file/directory; required without --input-manifest.",
    )
    parser.add_argument(
        "--json-input",
        "--metadata-input",
        dest="json_input",
        default=None,
        help="Skim-report JSON file/directory; required without --input-manifest.",
    )
    parser.add_argument(
        "--additional-metadata-input",
        "--extra-metadata-input",
        dest="additional_metadata_inputs",
        action="append",
        default=[],
        help="Extra metadata JSON; repeatable, with later inputs taking precedence.",
    )
    parser.add_argument("--input-files-file", help="Text file listing ROOT files, one per line.")
    parser.add_argument("--input-manifest", help="Validation manifest with known-good ROOT/JSON files.")
    parser.add_argument("--dataset-name", "--dataset", dest="dataset_name", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument(
        "--systematics",
        nargs="+",
        default=["Central"],
        help="Systematic keys/groups, e.g. Central, JERC Muon PU, or all.",
    )
    parser.add_argument("--list-systematics", action="store_true", help="List systematics and exit.")
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Process the first N valid ROOT files; normalization remains dataset-wide.",
    )
    parser.add_argument(
        "--file-shard",
        default=None,
        help=(
            "Process only shard INDEX of TOTAL, as 'INDEX/TOTAL' with "
            "1 <= INDEX <= TOTAL. Files are distributed round-robin, so the "
            "shards are size-balanced. The normalization denominator stays "
            "dataset-wide, hence hadd-ing every shard reproduces the "
            "unsharded output exactly."
        ),
    )
    parser.add_argument(
        "--input-file-batch-size",
        type=int,
        default=None,
        help=(
            "Maximum input ROOT files materialized together. Each batch is "
            "merged into the final output. By default DNN jobs use one file "
            "per batch and non-DNN jobs process all files together."
        ),
    )
    parser.add_argument("--rdf-threads", type=int, default=1, help="RDataFrame worker threads.")
    parser.add_argument(
        "--systematic-batch-size",
        type=int,
        default=2,
        help="Maximum systematic variations per RDF event pass.",
    )
    parser.add_argument(
        "--variable-batch-size",
        type=int,
        default=5,
        help=(
            "Maximum variables booked in one RDF event pass. Smaller batches "
            "reduce peak memory at the cost of additional event passes."
        ),
    )
    parser.add_argument("--variables", nargs="+")
    parser.add_argument("--no-region-sample-routing", dest="region_sample_routing", action="store_false", help="Explicitly allow DY/EWK outside their default generated mass region")
    parser.add_argument("--mass-regions", nargs="+", default=["mass_inclusive", "Z_sideband", "Signal_Fit"])
    parser.add_argument("--categories", nargs="+", default=["baseline", "ggF", "VBF"])
    parser.add_argument("--additional-cuts", default=None)
    parser.add_argument(
        "--disable-jet-horn-veto",
        action="store_true",
        help="Disable the era jet-horn veto without changing selections.yaml.",
    )
    parser.add_argument("--dryrun", action="store_true")
    parser.add_argument(
        "--dy-jet-components",
        "--jet-gen-components",
        "--pu-hard-jet-components",
        dest="dy_jet_components",
        action="store_true",
        help="Write exclusive reco/gen-matched 0J, 1J and 2J component files.",
    )
    parser.add_argument(
        "--jet-gen-component-processes",
        "--pu-hard-processes",
        nargs="+",
        default=["DY", "EWK"],
        help="Processes eligible for jet-component splitting.",
    )
    parser.add_argument(
        "--all-mc-jet-components", action="store_true",
        help="Split every MC process by reco multiplicity and gen matching.",
    )
    parser.add_argument(
        "--vbf-eta-regions",
        "--eta-components",
        action="store_true",
        help="Split VBF and its jet components into incl/CC/CF/FF eta regions.",
    )
    parser.add_argument(
        "--vbf-pair-definition",
        choices=["maxmjj", "hardest", "leading"],
        default="maxmjj",
        help="VBF jet pair: maxmjj (production), hardest (hardest preselected pair "
             "passing the cuts) or leading (two leading jets, both must pass the preselection); "
             "see common/vbf_pair_definition.py.",
    )
    parser.add_argument(
        "--dnn-model-set",
        choices=["updated", "legacy"],
        default="updated",
        help="DNN payload generation to use.",
    )
    parser.add_argument("--shift-z-sideband-dnn-mass", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    workflow_manifest = None
    if args.input_manifest:
        if not os.path.isfile(args.input_manifest):
            raise FileNotFoundError(f"Input manifest does not exist: {args.input_manifest}")
        workflow_manifest = read_manifest(args.input_manifest)
        if workflow_manifest.get("era") != args.era or workflow_manifest.get("dataset") != args.dataset_name:
            raise ValueError(
                "Input manifest era/dataset does not match histogram job: "
                f"manifest=({workflow_manifest.get('era')}, "
                f"{workflow_manifest.get('dataset')}), "
                f"requested=({args.era}, {args.dataset_name})"
            )
        manifest_stage = workflow_manifest.get("stage")
        if manifest_stage == "validation":
            if "valid_root_files" not in workflow_manifest:
                raise ValueError("Validation manifest is missing 'valid_root_files'")
            if "valid_json_files" not in workflow_manifest:
                raise ValueError("Validation manifest is missing 'valid_json_files'")
            manifest_has_no_inputs = (
                not workflow_manifest.get("valid_root_files", [])
                and not workflow_manifest.get("valid_json_files", [])
                and not workflow_manifest.get("invalid_root_files", [])
                and not workflow_manifest.get("invalid_json_files", [])
            )
            empty_input_failure = (
                manifest_has_no_inputs
                and workflow_manifest.get("failures")
                == ["no ROOT files or normalization JSON reports were discovered"]
            )
            if (
                workflow_manifest.get("status", "passed") != "passed"
                and not empty_input_failure
            ):
                invalid_roots = len(workflow_manifest.get("invalid_root_files", []))
                valid_roots = len(workflow_manifest.get("valid_root_files", []))
                raise RuntimeError(
                    "Refusing failed validation manifest "
                    f"{args.input_manifest}: {valid_roots} valid ROOT file(s), "
                    f"{invalid_roots} invalid ROOT file(s). Rerun validation "
                    "after completing the skim production."
                )
            if workflow_manifest.get("status", "passed") != "passed":
                print(
                    "[WARNING] Validation failed because no ROOT or JSON inputs "
                    "were discovered: producing an empty histogram output."
                )
            # The manifest is the source of truth. Folder arguments, when also
            # supplied, are intentionally ignored for the validated file lists.
            args.root_input = workflow_manifest.get("root_input")
            args.json_input = workflow_manifest.get("json_input")
            args.metadata_inputs = unique_metadata_inputs([
                *workflow_manifest["valid_json_files"], *args.additional_metadata_inputs
            ])
        else:
            raise ValueError(f"Unsupported histogram input manifest stage: {manifest_stage}")
    else:
        if not args.root_input:
            parser.error("--root-input is required when --input-manifest is not provided")
        if not args.json_input:
            parser.error("--json-input is required when --input-manifest is not provided")
        # Standalone mode: consume exactly the ROOT and JSON inputs supplied by
        # the caller. Validation is a separate workflow stage.
        args.metadata_inputs = expand_metadata_inputs([args.json_input, *args.additional_metadata_inputs])
    start_time = time.time()
    if args.max_files is not None and args.max_files < 1:
        raise ValueError("--max-files must be >= 1")
    if args.file_shard is not None:
        args.file_shard = parse_file_shard(args.file_shard)
    if args.rdf_threads < 1:
        raise ValueError("--rdf-threads must be >= 1")
    if args.systematic_batch_size < 1:
        raise ValueError("--systematic-batch-size must be >= 1")
    if args.variable_batch_size < 1:
        raise ValueError("--variable-batch-size must be >= 1")
    if args.input_file_batch_size is not None and args.input_file_batch_size < 1:
        raise ValueError("--input-file-batch-size must be >= 1")
    if args.rdf_threads > 1:
        ROOT.EnableImplicitMT(args.rdf_threads)
        print(f"[INFO] Enabled ROOT implicit multithreading with {args.rdf_threads} threads")
    analysis_path = os.environ["ANALYSIS_PATH"]
    cfg_dir = os.path.join(analysis_path, "config", args.era)
    main_cfg = utilities.get_config(os.path.join(cfg_dir, "maincfg.yaml"))
    samples_cfg = utilities.get_config(os.path.join(cfg_dir, "samples.yaml"))
    dataset_cfg = samples_cfg.get(args.dataset_name, {})
    is_data = dataset_cfg.get("is_data", False) or "data" in args.dataset_name.lower()
    # Validation is the authority for usable inputs. Every ROOT reaching this
    # stage must succeed, while all valid JSON reports in the manifest define
    # the full MC normalization denominator.
    requested_systematics = parse_requested_systematics(args.systematics)
    request_all = any(name.lower() == "all" for name in requested_systematics)
    request_central_only = {name.lower() for name in requested_systematics} <= {"central", "nominal"}
    systematics_mode = "central" if request_central_only else "all"
    if is_data and systematics_mode != "central":
        print(
            f"[SKIP] Dataset {args.dataset_name} is data: non-central "
            "systematic outputs are not produced."
        )
        sys.exit(0)
    empty_validated_input = bool(args.input_manifest and not workflow_manifest.get("valid_root_files", []))
    if empty_validated_input:
        # A failed/empty validation manifest has no events or variation
        # metadata from which systematic templates can be built.  Still emit
        # one nominal empty histogram so downstream campaign bookkeeping has a
        # valid ROOT output instead of failing while resolving (for example)
        # the QCDScale alias.
        print(
            "[WARNING] Validation manifest contains no valid ROOT files: "
            "producing Central empty histograms only."
        )
        requested_systematics = ["Central"]
        request_all = False
        request_central_only = True
        systematics_mode = "central"
    if args.list_systematics:
        request_all = True
        systematics_mode = "all"
    print(systematics_mode)
    args.systematics_mode = systematics_mode
    process_cfg = utilities.get_config(os.path.join(cfg_dir, "process_names.yaml"))
    args.process_name = utilities.process_from_dataset(process_cfg, args.dataset_name) or args.dataset_name
    if args.dy_jet_components and args.process_name.endswith("_nonStitched"):
        stitched_process = args.process_name.removesuffix("_nonStitched")
        stitched_entry = process_cfg.get(stitched_process, {})
        stitched_datasets = [
            *(stitched_entry.get("datasets", []) or []),
            *(stitched_entry.get("sub_processes", []) or []),
        ]
        if args.dataset_name in stitched_datasets:
            print(
                f"[INFO] Jet components: using stitched process "
                f"{stitched_process} instead of {args.process_name}."
            )
            args.process_name = stitched_process
    process_entry = process_cfg.get(args.process_name, {})
    sel_cfg = utilities.get_config(os.path.join(cfg_dir, "selections.yaml"))
    if args.disable_jet_horn_veto:
        sel_cfg["jet_horn_veto_expr"] = "(abs(v_ops::eta(Jet_p4)) < 0)"
    syst_cfg = utilities.get_config(os.path.join(cfg_dir, "systematics.yaml"))
    hist_cfg = utilities.get_config(os.path.join(analysis_path, "config", "plot", "histograms.yaml"))
    if args.dy_jet_components:
        if is_data:
            print(
                f"[INFO] Dataset {args.dataset_name} is data: producing "
                "reconstructed-jet categories without generator matching."
            )
            args.dy_jet_components = False
        elif args.all_mc_jet_components:
            print(f"[INFO] Splitting MC dataset {args.dataset_name} into jet components.")
        elif "split_jet_components" in process_entry:
            args.dy_jet_components = bool(process_entry["split_jet_components"])
            print(
                f"[INFO] Dataset {args.dataset_name} (process {args.process_name}): "
                f"split_jet_components={args.dy_jet_components} from process_names.yaml."
            )
        elif not jet_components_enabled_for_dataset(
                args.jet_gen_component_processes,
                args.dataset_name,
                args.process_name,
                is_signal=bool(process_entry.get("is_signal", False)),
            ):
            print(
                f"[INFO] Dataset {args.dataset_name} (process {args.process_name}) "
                "is outside --jet-gen-component-processes: producing normal "
                "histograms."
            )
            args.dy_jet_components = False
    if args.dy_jet_components:
        args.pu_hard_requested_categories = tuple(args.categories)
        sel_cfg = add_jet_component_categories(
            sel_cfg,
            include_vbf_eta_regions=args.vbf_eta_regions,
            requested_categories=args.pu_hard_requested_categories,
        )
        args.categories = list(
            expanded_jet_component_categories(
                include_vbf_eta_regions=args.vbf_eta_regions,
                requested_categories=args.pu_hard_requested_categories,
            )
        )
        requested_variables = list(args.variables if args.variables is not None else main_cfg["variables"])
        args.vbf_component_variables = requested_variables.copy()
        vars_to_add = ["eta_signed_vs_pt_leadingjet", "eta_signed_vs_pt_subleadingjet"]
        args.variables = list(dict.fromkeys([*requested_variables, *vars_to_add]))
    else:
        args.vbf_component_variables = []
        # After the component decision: data and the MC processes without jet
        # components must get the eta categories too, also when
        # --pu-hard-jet-components was requested for the campaign.
        if args.vbf_eta_regions:
            sel_cfg = add_vbf_eta_region_categories(sel_cfg)
            args.categories = [f"VBF_eta_{region}" for region in VBF_ETA_REGIONS]
    masses_regions = sel_cfg["masses_regions"]
    categories = sel_cfg["categories"]
    masses_regions_list = args.mass_regions
    categories_list = args.categories
    vars_to_make_hist = list(dict.fromkeys(args.variables or main_cfg["variables"]))
    dnn_payloads = sorted(
        {
            variable.rsplit("_NNOutput", 1)[0]
            for variable in vars_to_make_hist
            if variable.endswith("_NNOutput")
        }
    )
    btag_algo = main_cfg.get("bTagAlgo", "PNet")
    syst_cfg = configure_available_qcd_scale(syst_cfg, args.metadata_inputs, is_data, systematics_mode)
    systs_to_run = get_systs_to_run(syst_cfg, systematics_mode)
    if not request_all:
        requested_systematics = ["Central" if name.lower() == "nominal" else name for name in requested_systematics]
        systs_to_run = filter_systs_to_run(systs_to_run, requested_systematics)
    validate_systematic_isolation(systs_to_run)
    # Selection construction used to receive every object systematic from the
    # YAML even after the requested set had been filtered.  Asking for JERC
    # could therefore compile MuonScale expressions such as
    # m_mumu_FSR_scale_up.  Keep only the object-systematic families that are
    # actually present in systs_to_run; weight-only variations remain central
    # selections as intended.
    active_object_systematics = {"Central"}
    for base_name in syst_cfg.get("systematics", {}):
        if base_name == "Central":
            continue
        if any(
            run_name in systs_to_run
            for run_name in (f"{base_name}Up", f"{base_name}Down")
        ) or (
            base_name == "JER"
            and any(info.get("jer_component") for info in systs_to_run.values())
        ):
            active_object_systematics.add(base_name)
    syst_cfg = copy.deepcopy(syst_cfg)
    syst_cfg["systematics"] = {
        name: info
        for name, info in syst_cfg.get("systematics", {}).items()
        if name in active_object_systematics
    }
    active_jer_components = {
        info["jer_component"]
        for info in systs_to_run.values()
        if info.get("jer_component")
    }
    if active_jer_components and "JER" in syst_cfg["systematics"]:
        jer_info = copy.deepcopy(syst_cfg["systematics"]["JER"])
        jer_info["components"] = [
            component
            for component in jer_info.get("components", ())
            if component in active_jer_components
        ]
        syst_cfg["systematics"]["JER"] = jer_info
    if args.list_systematics:
        for syst_name in systs_to_run:
            print(syst_name)
        sys.exit(0)
    if empty_validated_input:
        dataset_seg_dict = {}
        print(
            "[INFO] No validated inputs: skipping normalization metadata "
            "loading for empty histogram production."
        )
    elif is_data:
        dataset_seg_dict = {}
        print(f"[INFO] Dataset {args.dataset_name} is data: skipping segmentation metadata loading.")
    else:
        metadata_started = time.perf_counter()
        print(f"[INFO] Loading normalization metadata from {len(args.metadata_inputs)} JSON input(s)...")
        dataset_seg_dict = get_combined_segmentation_dict(args.metadata_inputs)
        if not dataset_seg_dict:
            raise RuntimeError(
                "No normalization denominator could be read from --json-input. "
                "Use the report JSON paired with the selected skim ROOT file, "
                "or a directory containing valid reports."
            )
        profile_log(None, "central dataset metadata loading", metadata_started)
    dataset_qcd_scale_seg_dicts = {}
    if not is_data and systematics_mode != "central" and syst_cfg.get("qcd_scale", {}).get("enabled", False):
        qcd_metadata_started = time.perf_counter()
        for point in get_qcd_scale_points(syst_cfg["qcd_scale"]):
            point_name = point["name"]
            dataset_qcd_scale_seg_dicts[point_name] = get_qcd_scale_segmentation_dict(
                args.metadata_inputs, point_name, fallback_to_initial=False
            )
        profile_log(None, "QCD-scale dataset metadata loading", qcd_metadata_started)
    print(
        f"[INFO] Dataset metadata ready: {len(dataset_seg_dict)} central "
        f"entries and {len(dataset_qcd_scale_seg_dicts)} QCD-scale dictionaries"
    )
    if args.input_manifest:
        all_root_files = workflow_manifest["valid_root_files"]
    elif args.input_files_file:
        with open(args.input_files_file) as input_files_handle:
            all_root_files = [
                line.strip()
                for line in input_files_handle
                if line.strip() and not line.lstrip().startswith("#")
            ]
    else:
        all_root_files = list_root_files(args.root_input)
    if not args.input_manifest:
        validation_results = [validate_file((path, "Events")) for path in all_root_files]
        invalid_inputs = [(path, reason) for path, is_valid, reason in validation_results if not is_valid]
        if invalid_inputs:
            details = "\n".join(f"  {path}: {reason}" for path, reason in invalid_inputs)
            raise RuntimeError(f"Invalid histogram input ROOT file(s):\n{details}")
        all_root_files = [path for path, _, _ in validation_results]
    # Validation is external. Files from a manifest are already known-good;
    # standalone inputs are consumed exactly as supplied/discovered.
    valid_root_files = [os.path.abspath(path) for path in all_root_files]
    if args.max_files is not None:
        original_file_count = len(valid_root_files)
        valid_root_files = valid_root_files[: args.max_files]
        print(f"[TEST MODE] Processing {len(valid_root_files)}/{original_file_count} valid ROOT files.")
    if args.file_shard is not None:
        shard_index, shard_total = args.file_shard
        sharded_file_count = len(valid_root_files)
        valid_root_files = valid_root_files[shard_index - 1 :: shard_total]
        print(
            f"[SHARD {shard_index}/{shard_total}] Processing "
            f"{len(valid_root_files)}/{sharded_file_count} valid ROOT files."
        )
    if not valid_root_files:
        print("[WARNING] No validated ROOT files. Producing empty histograms.")
    if args.dryrun:
        print(f"[DRYRUN] Input files: {len(valid_root_files)}")
        print(f"[DRYRUN] Segmentation entries: {len(dataset_seg_dict)}")
        print(f"[DRYRUN] QCD-scale dictionaries: {len(dataset_qcd_scale_seg_dicts)}")
        for input_file in valid_root_files:
            print(f"  {input_file}")
        print("\n[DRYRUN] Exiting.")
        sys.exit(0)
    safe_mkdir(os.path.dirname(args.output_file))
    if os.path.exists(args.output_file):
        print(f"[INFO] Removing existing output file: {args.output_file}")
        os.remove(args.output_file)
    input_file_batch_size = args.input_file_batch_size or (
        1 if dnn_payloads else max(1, len(valid_root_files))
    )
    input_batches = (
        batch_list(valid_root_files, input_file_batch_size)
        if valid_root_files
        else [[]]
    )
    final_output_file = args.output_file
    batch_output_files = []
    # The per-batch parts and their merge stay on local disk (the job scratch on
    # Condor): hundreds of small writes and re-reads through EOS FUSE are slow, and
    # a single unreadable part used to lose a whole multi-hour job.  Only the merged
    # file is copied to the output location.
    scratch = None
    if len(input_batches) > 1:
        import tempfile
        scratch = tempfile.mkdtemp(
            prefix="hist_parts_",
            dir=os.environ.get("_CONDOR_SCRATCH_DIR") or os.environ.get("TMPDIR"),
        )
    try:
        for batch_index, input_batch in enumerate(input_batches, start=1):
            if len(input_batches) == 1:
                batch_output = final_output_file
            else:
                batch_output = os.path.join(
                    scratch,
                    f"{os.path.basename(final_output_file)}.part{batch_index:04d}.root",
                )
            args.output_file = batch_output
            batch_output_files.append(batch_output)
            print(
                f"[INPUT BATCH {batch_index}/{len(input_batches)}] "
                f"Processing {len(input_batch)} ROOT file(s)"
            )
            produce_histograms((
                input_batch,
                dataset_seg_dict,
                dataset_qcd_scale_seg_dicts,
                args,
                is_data,
                sel_cfg,
                syst_cfg,
                vars_to_make_hist,
                masses_regions,
                masses_regions_list,
                categories,
                categories_list,
                hist_cfg,
                systs_to_run,
                dnn_payloads,
                btag_algo,
            ))
        args.output_file = final_output_file
        if len(batch_output_files) > 1:
            merged_local = os.path.join(scratch, os.path.basename(final_output_file))
            merger = ROOT.TFileMerger(False, False)
            merger.OutputFile(merged_local, "RECREATE")
            for batch_output in batch_output_files:
                if not merger.AddFile(batch_output):
                    raise RuntimeError(
                        f"Could not add input-batch output: {batch_output}"
                    )
            if not merger.Merge():
                raise RuntimeError(
                    f"Could not merge input batches into {merged_local}"
                )
            shutil.copyfile(merged_local, final_output_file)
            if os.path.getsize(final_output_file) != os.path.getsize(merged_local):
                raise RuntimeError(f"Incomplete copy of the merged output: {final_output_file}")
    finally:
        args.output_file = final_output_file
        if scratch is not None:
            shutil.rmtree(scratch, ignore_errors=True)
    # An EOS write error can leave a truncated file that ROOT silently "recovers" when
    # reopened: fail instead, so that Condor retries and the check sees it missing.
    assert_complete_root_file(args.output_file)
    output = ROOT.TFile.Open(args.output_file, "UPDATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not reopen output file: {args.output_file}")
    if systematics_mode != "central":
        write_qcd_scale_variations(
            output,
            syst_cfg,
            vars_to_make_hist,
            masses_regions_list,
            categories_list,
            args.era,
            args.process_name,
        )
    output.Close()
    # The UPDATE above writes to EOS again: check before the DY split reads it back.
    assert_complete_root_file(args.output_file)
    written_files = [args.output_file]
    if args.dy_jet_components:
        written_files += split_dy_jet_component_outputs(
            args.output_file,
            masses_regions_list,
            process_label=args.process_name,
            include_vbf_eta_regions=args.vbf_eta_regions,
            requested_categories=args.pu_hard_requested_categories,
        )
    for written in written_files:
        assert_complete_root_file(str(written))
    execution_time = time.time() - start_time
    print("\n" + "=" * 80)
    print("[INFO] Histogram production completed successfully.")
    print(f"[INFO] Output file: {args.output_file}")
    print(f"[INFO] Input files: {len(valid_root_files)}")
    print(f"[INFO] Execution time:    {execution_time:.2f} s")
    print("=" * 80 + "\n")


if __name__ == "__main__":

    main()
