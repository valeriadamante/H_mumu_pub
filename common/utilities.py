"""Shared helpers for skim, validation, histograms, plotting, and campaigns.

Sections: configuration/datasets, files/reports/chunking, validation/manifests,
ROOT runtime, and histogram binning. ROOT is imported only by helpers using it,
so campaign and submission commands also work outside the CMSSW environment.
Dataframe physics preparation lives in common.prepare_rdf.
"""

import array
import ast
import bisect
from collections import OrderedDict
from datetime import datetime, timezone
from enum import Enum
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import yaml

_ROOT_RUNTIME_INITIALIZED = False
rootAnaPathSet = False
SCHEMA_VERSION = 1


# Configuration and dataset selection
# -----------------------------------

class WorkingPointsbTag(Enum):
    Loose = 1
    Medium = 2
    Tight = 3


def generate_enum_class(cls):
    enum_string = "enum class {} : int {{\n".format(cls.__name__)
    for item in cls:
        enum_string += "    {} = {},\n".format(item.name, item.value)
    enum_string += "};"
    return enum_string


def _resolve_config_path(yaml_file):
    path = Path(yaml_file)
    if path.exists() or path.is_absolute():
        return path
    return Path(__file__).resolve().parents[1] / path


def get_config(yaml_file):
    yaml_path = _resolve_config_path(yaml_file)
    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def process_from_dataset(process_cfg, dataset_name):
    """
    Iterates over process_names.yaml to find which process a dataset name belongs to.
    """
    for process, entry in process_cfg.items():
        if not isinstance(entry, dict):
            continue
        # A process with all datasets commented out leaves "datasets:" as None,
        # not as an empty list: without this guard .extend blew up with
        # AttributeError and stopped the producer on any dataset coming later
        # in the file. And it is a new list, not the one of the
        # configuration: extend() on the original glued the sub_processes
        # into it, polluting the shared config at every call.
        datasets_list = list(entry.get("datasets") or [])
        datasets_list += list(entry.get("sub_processes") or [])
        if dataset_name in datasets_list:
            return process
    return None


def GetObservablesCols(obs_name, is_data, nano_version="v12"):
    col_to_save_path = _resolve_config_path("config/col_to_save.yaml")
    with open(col_to_save_path, "r") as f:
        config_text = f.read()
    columns_literal = config_text.split("\ndef GetObservablesCols", 1)[0].strip()
    observables = ast.literal_eval(columns_literal)
    obs_to_store = []
    if obs_name not in observables.keys():
        raise RuntimeError(f"Invalid observable name, not found in keys {obs_name}")
    obs_dict = observables[obs_name]
    if "base" not in obs_dict.keys():
        raise RuntimeError(f"base key not found in {obs_dict}")
    obs_to_store.extend(obs_dict["base"])
    if nano_version in obs_dict.keys():
        obs_to_store.extend(obs_dict[nano_version])
    if not is_data and "MC" in obs_dict.keys():
        obs_to_store.extend(obs_dict["MC"])
    if "additional" in obs_dict.keys():
        obs_to_store.extend(obs_dict["additional"])
    return obs_to_store


def _load_yaml(path):
    with open(path) as handle:
        return yaml.safe_load(handle) or {}


def resolve_dataset_selection(analysis_path, era):
    """Resolve the canonical dataset/process selection from skim_cfg.yaml."""
    config_dir = Path(analysis_path) / "config" / era
    skim_cfg = _load_yaml(config_dir / "skim_cfg.yaml")
    processes_cfg = _load_yaml(config_dir / "process_names.yaml")
    samples_cfg = _load_yaml(config_dir / "samples.yaml")
    # samples_withfiles.yaml is generated (hmumu skim files) and not in git; it
    # only feeds the informational datasets_missing_filelist, so it is optional
    # for everything after the skims (hadd, checks, plots).
    files_path = config_dir / "samples_withfiles.yaml"
    samples_with_files = _load_yaml(files_path) if files_path.is_file() else None

    selected_processes = skim_cfg.get("process_to_select", []) or []
    whitelist = skim_cfg.get("datasets_whitelist", []) or []
    excluded = set(skim_cfg.get("datasets_exclude", []) or [])
    if not isinstance(selected_processes, list) or not isinstance(whitelist, list):
        raise ValueError(
            "skim_cfg.yaml process_to_select and datasets_whitelist must be lists"
        )

    process_datasets = OrderedDict()
    datasets = list(whitelist)
    for process in selected_processes:
        if process not in processes_cfg:
            raise KeyError(
                f"Process {process!r} from skim_cfg.yaml is missing from "
                f"config/{era}/process_names.yaml"
            )
        process_info = processes_cfg[process] or {}
        members = [
            *(process_info.get("datasets", []) or []),
            *(process_info.get("sub_processes", []) or []),
        ]
        members = list(dict.fromkeys(name for name in members if name not in excluded))
        process_datasets[process] = members
        datasets.extend(members)

    datasets = list(dict.fromkeys(name for name in datasets if name not in excluded))
    missing_samples = [name for name in datasets if name not in samples_cfg]
    if missing_samples:
        raise KeyError(
            f"Selected dataset(s) missing from config/{era}/samples.yaml: "
            + ", ".join(missing_samples)
        )

    return {
        "era": era,
        "datasets": datasets,
        "processes": list(process_datasets),
        "process_datasets": process_datasets,
        "datasets_whitelist": whitelist,
        "datasets_exclude": sorted(excluded),
        "datasets_missing_filelist": [] if samples_with_files is None else [
            name
            for name in datasets
            if not isinstance(samples_with_files.get(name), dict)
            or not isinstance(samples_with_files[name].get("filelist"), list)
        ],
    }


def load_routing(path):
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def era_policy(config, era):
    for policy in config["mass_region_routing"].values():
        if era in policy.get("eras", []):
            return policy
    raise ValueError(f"No mass-region sample routing configured for {era}")


def groups_for_region(config, era, mass_region):
    policy = era_policy(config, era)
    key = mass_region if mass_region in ("Signal_Fit", "H_sideband") else "sidebands"
    return tuple(policy[key])


def separate_groups(config, era):
    return tuple(era_policy(config, era).get("separate", []))


def jet_gen_component_processes(config):
    return tuple(config["jet_gen_components"]["enabled_processes"])


def dataset_region_allowed(dataset, region):
    """Use the matching generated mass window for DY/EWK, including FlashSim."""
    name = dataset.lower().replace("_", "")
    if not name.startswith(("dy", "ewk")):
        return True
    restricted = "105to160" in name
    return restricted == (region in ("Signal_Fit", "H_sideband"))


def production_samples(analysis_path, era, group):
    """Canonical nominal samples and configured FlashSim counterparts."""
    if group == 'FlashSim':
        candidates = {name for subset in ('signals', 'region_higgs', 'region_inclusive', 'flash_backgrounds')
                      for name in production_samples(analysis_path, era, subset) if 'flashsim' in name.lower()}
        samples = _load_yaml(Path(analysis_path) / 'config' / era / 'samples.yaml')
        return sorted(name for name in candidates
                      if any('Sep2026' in path for path in samples.get(name, {}).get('nanoAOD', [])))
    folder = Path(analysis_path) / 'config' / era
    samples = _load_yaml(folder / 'samples.yaml')
    processes = _load_yaml(folder / 'process_names.yaml')
    skim = _load_yaml(folder / 'skim_cfg.yaml')
    excluded = set(skim.get('datasets_exclude', []) or [])
    selected = []
    for process, entry in processes.items():
        names = [*(entry.get('datasets', []) or []), *(entry.get('sub_processes', []) or [])]
        for name in names:
            lower = name.lower()
            signal = lower.startswith(('glugluh', 'vbfh')) and 'to2mu' in lower
            if group == 'signals':
                keep = signal and not any(x in lower for x in ('120', '130', 'tune', 'minnlo'))
            elif group in ('region_higgs', 'region_inclusive'):
                nominal = process in ('DY', 'DYto2Mu_MLL105To160', 'EWK',
                    'EWK_2Mu2J_MLL_105to160_herwig', 'EWK_2Mu2J_MLL_105to160_pythia',
                    'EWK_LLJJ_herwig', 'EWK_LLJJ_pythia')
                keep = name.lower().startswith(('dy', 'ewk')) and (nominal or 'flashsim' in lower)
                keep = keep and dataset_region_allowed(name, 'H_sideband' if group == 'region_higgs' else 'Z_sideband')
            elif group == 'flash_backgrounds':
                keep = 'flashsim' in lower and not signal and not lower.startswith(('dy', 'ewk'))
            else:
                raise ValueError(group)
            if keep and name in samples and name not in excluded:
                selected.append(name)
    return sorted(set(selected))


# Files, skim reports, and input chunks
# -------------------------------------

def list_root_files(path):
    if path.endswith(".root"):
        return [path]

    files = []
    for root, _, fnames in os.walk(path):
        for fname in fnames:
            if fname.endswith(".root"):
                files.append(os.path.join(root, fname))
    return sorted(files)


def extract_dataset_name(path):
    return os.path.basename(path.rstrip("/"))


def report_path_for_root(root_file):
    directory = os.path.dirname(root_file)
    stem = os.path.splitext(os.path.basename(root_file))[0]
    if stem.startswith("skim_") and stem[len("skim_"):].isdigit():
        return os.path.join(directory, f"report_{stem[len('skim_'):]}.json")
    return os.path.splitext(root_file)[0] + "_report.json"


def load_report_json(root_file):
    json_file = report_path_for_root(root_file)
    if not os.path.exists(json_file):
        print(f"[WARNING] Missing json: {json_file}")
        return {}
    with open(json_file) as f:
        return json.load(f)


def SaveReport(rdf, report, verbose=0):

    cuts = [c for c in report]
    report_json = {}
    if len(cuts) > 0:
        # iniziale
        initial = cuts[0].GetAll()
        rdf = rdf.Define("Report_Initial", f"{initial}")
        report_json["Initial"] = initial
        for c_id, cut in enumerate(cuts):
            cut_name = '_'.join(str(cut.GetName()).split(" "))
            passed = cut.GetPass()
            eff = cut.GetEff()
            report_json[cut_name] = {
                "pass": passed,
                # "eff": eff
            }
            if verbose > 0:
                print(
                    f"for the cut {cut.GetName()} "
                    f"there are {passed} events passed over {initial}, "
                    f"resulting in an efficiency of {eff}"
                )

    return rdf, report_json


def chunk_files_by_size(file_entries, target_bytes, max_files):
    """Create stable sequential chunks bounded by input bytes and file count."""
    if target_bytes <= 0:
        raise ValueError("target_bytes must be greater than zero")
    if max_files <= 0:
        raise ValueError("max_files must be greater than zero")

    chunks = []
    current = []
    current_size = 0
    for raw_entry in file_entries:
        entry = (
            {"path": raw_entry, "size": 0}
            if isinstance(raw_entry, str)
            else dict(raw_entry)
        )
        if "path" not in entry:
            raise ValueError(f"file entry has no path: {raw_entry!r}")
        entry_size = max(0, int(entry.get("size", 0)))
        if current and (
            current_size + entry_size > target_bytes or len(current) >= max_files
        ):
            chunks.append(current)
            current = []
            current_size = 0
        current.append(entry)
        current_size += entry_size
    if current:
        chunks.append(current)
    return chunks


def input_file_id(path):
    stem = Path(path).stem
    match = re.fullmatch(r"(?:skim|report)_(\d+)", stem)
    return match.group(1) if match else stem.removesuffix("_report")


def metadata_excluding_root_files(metadata_inputs, excluded_root_files):
    """Return JSON inputs excluding reports paired to rejected MC ROOT files."""
    excluded_ids = {input_file_id(path) for path in excluded_root_files}
    selected = []
    seen = set()
    for raw_path in metadata_inputs:
        path = Path(raw_path)
        candidates = path.rglob("*.json") if path.is_dir() else (path,)
        for candidate in candidates:
            normalized = str(candidate.resolve())
            if (
                candidate.suffix == ".json"
                and normalized not in seen
                and input_file_id(candidate) not in excluded_ids
            ):
                selected.append(normalized)
                seen.add(normalized)
    return selected


def normalize_stem(path):
    base = os.path.basename(str(path))
    for ext in (".root", ".json"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


def filter_seg_dict_for_files(seg_dict, root_files):
    valid_stems = {normalize_stem(f) for f in root_files}

    out = {}
    for k, v in seg_dict.items():
        if normalize_stem(k) in valid_stems:
            out[k] = v

    return out


def get_segmentation_dict(
    json_paths,
    node="gen",
    fallback_to_initial=True,
    warn_if_missing=True,
):
    global_segmentation = {}

    # A single path is a common caller input.  Treat it as one path rather
    # than iterating over the individual characters of the string.
    if isinstance(json_paths, (str, os.PathLike)):
        json_paths = [json_paths]

    for json_path in json_paths:
        # Skipping a report would drop its events from the sum of weights and
        # over-normalise the dataset, so an unreadable report is fatal.  EOS I/O
        # errors are often transient: retry before giving up.
        for attempt in range(5):
            try:
                with open(json_path) as json_file:
                    info = json.load(json_file)
                break
            except (OSError, json.JSONDecodeError) as error:
                print(f"[WARNING] Could not parse JSON file {json_path} "
                      f"(attempt {attempt + 1}/5): {error}")
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError(f"Unreadable normalization JSON {json_path}")

        node_dict = info.get(node)

        if isinstance(node_dict, dict):
            segmented_keys = [
                key
                for key, value in node_dict.items()
                if key != "total" and isinstance(value, dict)
            ]

            if segmented_keys:
                # Dataset with distinct denominators per segmentation.
                # The total is ignored because it is cumulative.
                for sub_key in segmented_keys:
                    sub_info = node_dict[sub_key]

                    selection = sub_info.get("selection")
                    value = sub_info.get("value", 0.0)

                    if not selection:
                        continue

                    try:
                        value = float(value)
                    except (TypeError, ValueError):
                        print(
                            f"[WARNING] Invalid value for node '{node}', "
                            f"segment '{sub_key}', file {json_path}: {value!r}"
                        )
                        continue

                    global_segmentation[selection] = (
                        global_segmentation.get(selection, 0.0) + value
                    )

            elif isinstance(node_dict.get("total"), dict):
                # Non-segmented dataset.
                total_info = node_dict["total"]

                selection = total_info.get("selection", "return true;")
                value = total_info.get("value", 0.0)

                try:
                    value = float(value)
                except (TypeError, ValueError):
                    print(
                        f"[WARNING] Invalid total value for node '{node}', "
                        f"file {json_path}: {value!r}"
                    )
                    continue

                global_segmentation[selection] = (
                    global_segmentation.get(selection, 0.0) + value
                )

            else:
                print(
                    f"[WARNING] Node '{node}' has no usable segmentation "
                    f"entries in {json_path}"
                )

        elif fallback_to_initial and "Initial" in info:
            initial = info["Initial"]

            try:
                if isinstance(initial, dict):
                    value = sum(float(item) for item in initial.values())
                else:
                    value = float(initial)
            except (TypeError, ValueError):
                print(
                    f"[WARNING] Invalid Initial value in "
                    f"{json_path}: {initial!r}"
                )
                continue

            global_segmentation["return true;"] = (
                global_segmentation.get("return true;", 0.0) + value
            )


    return global_segmentation


# Validation and manifests
# ------------------------

def discover_root_files(path):
    if path.endswith(".root"):
        return [os.path.abspath(path)]
    files = []
    for root, _, names in os.walk(path):
        files.extend(
            os.path.abspath(os.path.join(root, name))
            for name in names
            if name.endswith(".root")
        )
    return sorted(files)


def validate_file(task):
    """Validate one ROOT tree; suitable for multiprocessing workers."""
    path, tree_name, *retry_options = task
    retries = int(retry_options[0]) if retry_options else 1
    retry_delay = float(retry_options[1]) if len(retry_options) > 1 else 0.0
    import ROOT

    last_reason = "unknown validation error"
    for attempt in range(1, retries + 1):
        root_file = None
        try:
            root_file = ROOT.TFile.Open(path, "READ")
            if not root_file or root_file.IsZombie():
                last_reason = "cannot open file or zombie"
            else:
                tree = root_file.Get(tree_name)
                if not tree:
                    last_reason = f"missing tree '{tree_name}'"
                elif tree.GetEntries() == 0:
                    last_reason = f"empty tree '{tree_name}' (0 entries)"
                elif tree.GetListOfBranches().GetEntries() == 0:
                    last_reason = f"tree '{tree_name}' has no branches"
                else:
                    return path, True, ""
        except Exception as error:
            last_reason = repr(error)
        finally:
            if root_file:
                root_file.Close()
        if attempt < retries and retry_delay:
            time.sleep(retry_delay)
    return path, False, f"{last_reason} (failed after {retries} attempts)"


def atomic_write_lines(path, lines):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w") as handle:
        for line in lines:
            handle.write(f"{line}\n")
    os.replace(temporary, output)


def nonempty(path):
    candidate = Path(path)
    return candidate.is_file() and candidate.stat().st_size > 0


def validation_complete(path):
    if not nonempty(path):
        return False
    try:
        with open(path) as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return False
    return (
        manifest.get("stage") == "validation"
        and manifest.get("status") == "passed"
        and isinstance(manifest.get("valid_root_files"), list)
        and (
            bool(manifest["valid_root_files"])
            or bool(manifest.get("ignored_empty_root_files", []))
        )
    )


def histogram_complete(path):
    """Return true only for a readable ROOT file containing output objects."""
    if not nonempty(path):
        return False
    root_file = None
    try:
        import ROOT

        root_file = ROOT.TFile.Open(str(path), "READ")
        return bool(
            root_file
            and not root_file.IsZombie()
            and root_file.GetNkeys() > 0
        )
    except Exception:
        return False
    finally:
        if root_file:
            root_file.Close()


def stage_output_complete(stage, path):
    if stage not in {"validation", "histograms", "systematics"}:
        raise ValueError(f"Unknown workflow stage: {stage}")
    return (
        validation_complete(path)
        if stage == "validation"
        else histogram_complete(path)
    )


def read_manifest(path, expected_stage=None):
    with open(path) as handle:
        manifest = json.load(handle)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Unsupported manifest schema in {path}")
    if expected_stage and manifest.get("stage") != expected_stage:
        raise ValueError(
            f"Expected a {expected_stage!r} manifest, got "
            f"{manifest.get('stage')!r}: {path}"
        )
    return manifest


def write_manifest(path, stage, **payload):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "schema_version": SCHEMA_VERSION,
        "stage": stage,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **payload,
    }
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w") as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, output)
    return document


def ensure_validation_manifest(
    manifest_path,
    *,
    era,
    dataset,
    root_input=None,
    json_input=None,
    fallback_path,
    workers=8,
    retries=3,
    retry_delay=2.0,
):
    """Return a valid manifest, creating it from raw inputs when necessary."""
    candidate = Path(manifest_path) if manifest_path else Path(fallback_path)
    if candidate.is_file():
        manifest = read_manifest(candidate, "validation")
    else:
        if not root_input:
            raise ValueError(
                "No validation manifest found: --input is required for "
                "automatic validation"
            )
        json_input = json_input or root_input
        analysis_path = Path(os.environ["ANALYSIS_PATH"])
        command = [
            sys.executable,
            str(analysis_path / "analysis/validate_dataset.py"),
            "--era",
            era,
            "--dataset-name",
            dataset,
            "--root-input",
            root_input,
            "--json-input",
            json_input,
            "--output-manifest",
            str(candidate),
            "--workers",
            str(workers),
            "--retries",
            str(retries),
            "--retry-delay",
            str(retry_delay),
        ]
        print("[MANIFEST] Not found; running automatic validation", flush=True)
        print("[MANIFEST] " + " ".join(command), flush=True)
        subprocess.run(command, check=True, cwd=analysis_path)
        manifest = read_manifest(candidate, "validation")
    if manifest.get("status", "passed") != "passed":
        raise RuntimeError(f"Refusing failed validation manifest: {candidate}")
    if manifest["era"] != era or manifest["dataset"] != dataset:
        raise ValueError(f"Manifest era/dataset mismatch: {candidate}")
    return str(candidate.resolve()), manifest


def is_valid_root_file(filename, tree_name="Events"):
    import ROOT

    root_file = None
    try:
        root_file = ROOT.TFile.Open(filename)
        if not root_file or root_file.IsZombie():
            return False
        tree = root_file.Get(tree_name)
        if not tree:
            return False
        return tree.GetEntries() > 0
    except Exception:
        return False
    finally:
        if root_file:
            root_file.Close()


def get_valid_root_files(files, tree_name="Events"):
    valid_files = []
    for f in files:
        if is_valid_root_file(f, tree_name):
            valid_files.append(f)
        else:
            print(f"[WARNING] Skipping invalid file: {f}")
    return valid_files


# ROOT runtime and column helpers
# -------------------------------

def initialize_root_runtime(batch=True, thread_safe=True):
    """Initialize ROOT and the shared analysis header once per process."""
    import ROOT

    global _ROOT_RUNTIME_INITIALIZED
    if _ROOT_RUNTIME_INITIALIZED:
        return
    analysis_path = os.environ.setdefault(
        "ANALYSIS_PATH", str(Path(__file__).resolve().parents[1])
    )
    if batch:
        ROOT.gROOT.SetBatch(True)
    if thread_safe:
        ROOT.EnableThreadSafety()
    DeclareHeader(f"{analysis_path}/analysis/AnalysisTools.h")
    _ROOT_RUNTIME_INITIALIZED = True


def _column_names(df):
    return {str(col) for col in df.GetColumnNames()}


def _has_column(df, column):
    return column in _column_names(df)


def _define_if_missing(df, name, expression):
    if _has_column(df, name):
        return df
    return df.Define(name, expression)


def ListToVector(list, type="string"):
    import ROOT

    vec = ROOT.std.vector(type)()
    for item in list:
        vec.push_back(item)
    return vec


def mkdir_recursive(root_file, dir_path):
    if dir_path == "":
        return root_file
    current = root_file
    for folder in dir_path.split("/"):
        if not current.GetDirectory(folder):
            current.mkdir(folder)
        current = current.GetDirectory(folder)
    return current


def mkdir(file, path):
    import ROOT

    dir_names = path.split("/")
    current_dir = file
    for n, dir_name in enumerate(dir_names):
        dir_obj = current_dir.Get(dir_name)
        full_name = f"{file.GetPath()}" + "/".join(dir_names[:n])
        if dir_obj:
            if not dir_obj.IsA().InheritsFrom(ROOT.TDirectory.Class()):
                raise RuntimeError(
                    f"{dir_name} already exists in {full_name} and it is not a directory"
                )
        else:
            dir_obj = current_dir.mkdir(dir_name)
            if not dir_obj:

                raise RuntimeError(f"Failed to create {dir_name} in {full_name}")
        current_dir = dir_obj
    return current_dir


def DeclareHeader(header, verbose=0):
    import ROOT

    global rootAnaPathSet
    if not rootAnaPathSet:
        if verbose > 0:
            print(f'Adding "{os.environ["ANALYSIS_PATH"]}" to the ROOT include path')
        ROOT.gROOT.ProcessLine(".include " + os.environ["ANALYSIS_PATH"])
        rootAnaPathSet = True
    if verbose > 0:
        print(f'Including "{header}"')
    if not os.path.exists(header):
        raise RuntimeError(f'"{header}" does not exist')
    if not ROOT.gInterpreter.Declare(f'#include "{header}"'):
        raise RuntimeError(f"Failed to include {header}")
    if verbose > 0:
        print(f'Successfully included "{header}"')


# Histogram models, binning, and contents
# ---------------------------------------

def findBinEntry(hist_cfg_dict, var_name):
    """
    Match variable name against regex-based histogram config entries.
    """

    matches = []

    for pattern in hist_cfg_dict.keys():
        if re.fullmatch(pattern, var_name):
            matches.append(pattern)

    if not matches:
        raise KeyError(f"No histogram config pattern matches variable '{var_name}'")

    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous histogram config for '{var_name}': {matches}")

    return matches[0]


def findNewBins(hist_cfg_dict, var, **keys):
    cfg = hist_cfg_dict.get(var, {})
    if 'x_rebin' not in cfg: return cfg.get('x_bins', [])
    x_rebin = cfg['x_rebin']
    if isinstance(x_rebin, list): return x_rebin
    def recursive_search(d, remaining_keys):
        if isinstance(d, list): return d
        if not remaining_keys and isinstance(d, dict) and 'other' in d: return d['other']
        if not isinstance(d, dict): return None
        for k_name, k_value in remaining_keys.items():
            if k_value in d:
                found = recursive_search(d[k_value], {kk: vv for kk, vv in remaining_keys.items() if kk != k_name})
                if found is not None: return found
        return d.get('other') if isinstance(d, dict) else None
    return recursive_search(x_rebin, {k: v for k, v in keys.items() if v is not None}) or cfg.get('x_bins', [])


def GetModel(hist_cfg, var, dims, era=None):
    import ROOT

    THModel_Inputs = []
    var_entry = findBinEntry(hist_cfg, var)
    if dims == 1:
        variable_cfg = hist_cfg[var_entry]
        binning = variable_cfg.get("binning", {})
        x_bins = binning.get(era, variable_cfg["x_bins"])
        x_bins_vec = GetBinVec(x_bins)
        THModel_Inputs.append(x_bins_vec.size() - 1)
        THModel_Inputs.append(x_bins_vec.data())
        model = ROOT.RDF.TH1DModel("", "", *THModel_Inputs)
        # TH1DModel keeps the edge pointer rather than copying its storage.
        # Retain the vector for at least as long as the Python model wrapper.
        model._bin_vectors = [x_bins_vec]
        return model

    elif (dims == 2) or (dims == 3):
        list_var_bins_vec = []
        for var_nD in hist_cfg[var_entry]["var_list"]:
            var_bin_name = f"{var_nD}_bins"
            var_bins = (
                hist_cfg[var_entry][var_bin_name]
                if var_bin_name in hist_cfg[var_entry]
                else hist_cfg[var_nD]["x_bins"]
            )
            var_bins_vec = GetBinVec(var_bins)
            list_var_bins_vec.append(var_bins_vec)
            THModel_Inputs.append(var_bins_vec.size() - 1)
            THModel_Inputs.append(var_bins_vec.data())
        if dims == 2:
            model = ROOT.RDF.TH2DModel("", "", *THModel_Inputs)
            model._bin_vectors = list_var_bins_vec
            return model
        if dims == 3:
            model = ROOT.RDF.TH3DModel("", "", *THModel_Inputs)
            model._bin_vectors = list_var_bins_vec
            return model
            return model
    else:
        raise RuntimeError("nD histogram not implemented yet")
        # model = ROOT.RDF.THnDModel("", "", )

    return model


def getNewBins(bins):
    if isinstance(bins, list): return bins
    n_bins_str, bin_range = bins.split('|')
    start, stop = map(float, bin_range.split(':'))
    n_bins = int(n_bins_str)
    return [start + i * (stop - start) / n_bins for i in range(n_bins + 1)]


def GetBinVec(x_bins):
    import numpy as np

    if isinstance(x_bins, dict):
        return x_bins

    x_bins_vec = None
    if not isinstance(x_bins, list):
        n_bins, bin_range = x_bins.split("|")
        start, stop = bin_range.split(":")
        x_bins = np.linspace(float(start), float(stop), int(n_bins) + 1).tolist()
    x_bins_vec = ListToVector(x_bins, "float")
    return x_bins_vec


def AdaptBinningToHistogram(hist, desired_binning):
    axis = hist.GetXaxis()
    original_edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)]
    adapted = []
    for x in desired_binning:
        idx = bisect.bisect_left(original_edges, x)
        if idx == 0: closest = original_edges[0]
        elif idx == len(original_edges): closest = original_edges[-1]
        else:
            before, after = original_edges[idx - 1], original_edges[idx]
            closest = before if abs(x - before) < abs(x - after) else after
        adapted.append(closest)
    return sorted(set(adapted))


def FixNegativeContributions(histogram):
    orig_integral = histogram.Integral(0, histogram.GetNbinsX() + 1)
    if orig_integral < 0:
        print(f"Integral negative for {histogram.GetName()}")
        return False, "", ""
    for n in range(1, histogram.GetNbinsX() + 1):
        if histogram.GetBinContent(n) < 0:
            error = abs(histogram.GetBinContent(n))
            new_error = math.sqrt(error**2 + histogram.GetBinError(n)**2)
            histogram.SetBinContent(n, 0)
            histogram.SetBinError(n, new_error)
    if orig_integral > 0: histogram.Scale(1.0)
    return True, "", ""


def RebinHisto(hist_initial, new_binning, sample, wantOverflow=True, verbose=False):
    import ROOT

    adapted = AdaptBinningToHistogram(hist_initial, new_binning)
    if len(adapted) < 2:
        print(f"Rebinning histogram {hist_initial.GetName()} with {len(new_binning)-1} bins")
        print(hist_initial.GetName(), hist_initial.GetNbinsX(), hist_initial.GetXaxis().GetXmin(), hist_initial.GetXaxis().GetXmax())
        print("Initial binning:", [hist_initial.GetXaxis().GetBinLowEdge(i) for i in range(1, hist_initial.GetNbinsX() + 2)])
        print(adapted)
        raise RuntimeError("Adapted binning < 2 edges!")
    new_hist = hist_initial.Rebin(len(adapted) - 1, sample, array.array('d', adapted))
    if sample == 'data': new_hist.SetBinErrorOption(ROOT.TH1.kPoisson)
    if wantOverflow:
        n_final = new_hist.GetBinContent(new_hist.GetNbinsX())
        n_over = new_hist.GetBinContent(new_hist.GetNbinsX() + 1)
        new_hist.SetBinContent(new_hist.GetNbinsX(), n_final + n_over)
    FixNegativeContributions(new_hist)
    return new_hist


def is_valid_histogram(hist, check_overflow=True):
    """
    Check that the histogram is usable for plotting.
    """
    import numpy as np


    if hist is None:
        return False

    if not hist.InheritsFrom("TH1"):
        return False

    n_bins = hist.GetNbinsX()

    total_content = 0.0

    first_bin = 0 if check_overflow else 1
    last_bin  = n_bins + 1 if check_overflow else n_bins

    for ibin in range(first_bin, last_bin + 1):

        content = hist.GetBinContent(ibin)
        error   = hist.GetBinError(ibin)

        if not np.isfinite(content):
            return False

        if not np.isfinite(error):
            return False

        total_content += abs(content)

    return total_content > 0
