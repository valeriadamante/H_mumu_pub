import argparse
import copy
import os
import re
import shutil
import sys

import yaml

ANALYSIS_PATH = os.environ.get(
    "ANALYSIS_PATH",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
)
sys.path.append(ANALYSIS_PATH)

from common.utilities import resolve_dataset_selection
from common.jet_component_splitting import DY_COMPONENT_FILE_LABELS
from common.combine_processes import combine_datasets


def load_yaml_config(yaml_path):
    with open(yaml_path, "r") as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(f"Error parsing the YAML file: {exc}")
            return None


def hist_sum_value(histogram):
    total = histogram.sum(flow=True)
    return float(getattr(total, "value", total))


def set_hist_values(histogram, values):
    """Set bin contents for both weighted and unweighted hist storages."""
    view = histogram.view(flow=True)
    target = getattr(view, "value", view)
    target[...] = values


def dataset_file_candidates(input_dir, process, dataset, component_label=None):
    suffix_by_process = {
        "DYto2Mu_MLL105To160": ["_stitched", ""],
        "DYto2Mu_MLL105To160_nonStitched": ["_nonStitched", ""],
        "DYto2Mu_MLL105To160_VBFFiltered": ["_stitched", ""],
        "DYto2Mu_MLL105To160_VBFFiltered_nonStitched": ["_nonStitched", ""],
        "DYto2Mu_MLL105To160_FlashSim": ["_nonStitched", ""],
    }
    source_process = process
    suffixes = suffix_by_process.get(source_process, [""])
    # hist_maker prefixes split files with the concrete process name, except
    # for the canonical DY process where the historical DY_* labels are kept.
    # Examples:
    #   DY                         -> <dataset>_DY_VBF_Hard.root
    #   DYto2Mu_MLL105To160        -> <dataset>_DYto2Mu_MLL105To160_VBF_Hard.root
    candidates = []
    component_name = component_label.removeprefix("DY_") if component_label else None
    for suffix in suffixes:
        if component_name and source_process != "DY":
            # Current split files encode the concrete process variant inside
            # the component suffix, e.g.
            #   <dataset>_DYto2Mu_MLL105To160_nonStitched_0J.root
            candidates.append(
                os.path.join(
                    input_dir,
                    f"{dataset}_{source_process}{suffix}_{component_name}.root",
                )
            )
            # Retain compatibility with the older dataset-suffix layout.
            candidates.append(
                os.path.join(
                    input_dir,
                    f"{dataset}{suffix}_{source_process}_{component_name}.root",
                )
            )
        elif component_name:
            candidates.append(
                os.path.join(input_dir, f"{dataset}{suffix}_{component_label}.root")
            )
        else:
            candidates.append(os.path.join(input_dir, f"{dataset}{suffix}.root"))
    return list(dict.fromkeys(candidates))


def luminosity(era):
    """Integrated luminosity of an era (pb-1), from its maincfg.yaml."""
    return float(load_yaml_config(os.path.join("config", era, "maincfg.yaml"))["luminosity"])


def add_derived_systematics(era, output_dir):
    """Shape variations built from two samples (e.g. herwig vs pythia).

    Up/Down = target +- coefficient * scale * |alternative - nominal|.  Normally
    nominal, alternative and target are the same era and the target is the
    nominal sample.  ``target_process`` applies the difference of another pair
    (e.g. a private herwig/pythia pair to the official herwig sample).  With
    ``reference_era`` the difference is taken from that
    era (the sibling directory of output_dir) and applied to ``target_process``
    of this era, scaled by lumi(era)/lumi(reference) when ``scale: lumi``: for
    an era without the alternative sample.
    """
    import numpy as np
    import uproot

    systematics_cfg = load_yaml_config(
        os.path.join("config", era, "systematics.yaml")
    )
    derived_cfg = (systematics_cfg or {}).get("derived_systematics", {})

    for _, config in derived_cfg.items():
        nominal_process = config["nominal_process"]
        alternative_process = config["alternative_process"]
        from common.systematic_correlations import nuisance_name as correlated_nuisance_name
        nuisance_name = correlated_nuisance_name(config, era)
        coefficient = float(config.get("coefficient", 0.5))
        floor = float(config.get("floor", 0.0))
        reference_era = config.get("reference_era")
        reference_dir = output_dir
        scale = 1.0
        if reference_era:
            reference_dir = os.path.join(os.path.dirname(os.path.normpath(output_dir)), reference_era)
            scale = (luminosity(era) / luminosity(reference_era)
                     if config.get("scale") == "lumi" else float(config.get("scale", 1.0)))
        nominal_path = os.path.join(reference_dir, f"{nominal_process}.root")
        alternative_path = os.path.join(reference_dir, f"{alternative_process}.root")
        target_process = config.get("target_process", nominal_process)
        target_path = os.path.join(output_dir, f"{target_process}.root")

        missing = [path for path in dict.fromkeys((nominal_path, alternative_path, target_path))
                   if not os.path.exists(path)]
        if missing:
            print(f"[WARNING] Cannot build {nuisance_name}: missing {', '.join(missing)}")
            continue

        variations = {}
        with (
            uproot.open(nominal_path) as nominal_file,
            uproot.open(alternative_path) as alternative_file,
            uproot.open(target_path) as target_file,
        ):
            for key in target_file.keys(recursive=True):
                clean_key = key.split(";")[0]
                histogram_name = clean_key.rsplit("/", 1)[-1]
                if re.search(r"(Up|Down)$", histogram_name):
                    continue
                if clean_key not in nominal_file or clean_key not in alternative_file:
                    continue

                objects = (target_file[clean_key], nominal_file[clean_key], alternative_file[clean_key])
                if not all(hasattr(obj, "to_hist") for obj in objects):
                    continue

                target_hist = objects[0].to_hist()
                target_values = target_hist.values(flow=True)
                nominal_values = objects[1].to_hist().values(flow=True)
                alternative_values = objects[2].to_hist().values(flow=True)
                if not (target_values.shape == nominal_values.shape == alternative_values.shape):
                    raise RuntimeError(
                        f"Histogram shape mismatch for '{clean_key}' between "
                        f"{target_process}, {nominal_process} and {alternative_process}"
                    )

                half_difference = coefficient * scale * np.abs(
                    alternative_values - nominal_values
                )
                up_hist = copy.deepcopy(target_hist)
                down_hist = copy.deepcopy(target_hist)
                if np.any(target_values):
                    set_hist_values(up_hist, target_values + half_difference)
                    set_hist_values(
                        down_hist,
                        np.maximum(floor, target_values - half_difference),
                    )
                else:
                    # Keep empty nominal templates empty also for derived
                    # systematics.  The deep copies retain axes and metadata.
                    set_hist_values(up_hist, np.zeros_like(target_values))
                    set_hist_values(down_hist, np.zeros_like(target_values))

                directory = clean_key.rsplit("/", 1)[0] if "/" in clean_key else ""
                prefix = f"{directory}/" if directory else ""
                variations[
                    f"{prefix}{histogram_name}_{nuisance_name}Up"
                ] = up_hist
                variations[
                    f"{prefix}{histogram_name}_{nuisance_name}Down"
                ] = down_hist

        with uproot.update(target_path) as target_output:
            for key, histogram in variations.items():
                target_output[key] = histogram
        source = f" (from {reference_era}, x{scale:.4f})" if reference_era else ""
        print(
            f"   -> Added {len(variations)} {nuisance_name} templates "
            f"to {target_path}{source}"
        )


def write_histograms(output_file_path, histograms, attempts=3):
    """Write to a local file and then copy it: uproot.recreate directly on EOS
    FUSE does non-sequential writes and leaves unreadable histograms
    (R__unzip_header, 'unrecognized compression'), seen on 2026-09-24.
    After the copy, reread every key and retry if something is wrong.
    """
    import tempfile
    import uproot

    with tempfile.TemporaryDirectory(prefix="hadd_", dir=os.environ.get("TMPDIR")) as tmp:
        local = os.path.join(tmp, os.path.basename(output_file_path))
        with uproot.recreate(local) as output_root:
            for key, merged_hist in histograms.items():
                output_root[key] = merged_hist
        for attempt in range(1, attempts + 1):
            shutil.copyfile(local, output_file_path)
            try:
                with uproot.open(output_file_path) as check:
                    for key in histograms:
                        check[key].values()
                # uproot does not see everything: a file can open with uproot and
                # still be unreadable for ROOT (damaged header).
                import ROOT
                handle = ROOT.TFile.Open(output_file_path)
                if not handle or handle.IsZombie():
                    raise OSError("ROOT cannot open the file")
                handle.Close()
                return
            except Exception as error:
                print(f"   WARNING: {output_file_path} unreadable after copy "
                      f"(attempt {attempt}/{attempts}): {error}")
        raise RuntimeError(f"{output_file_path} is still unreadable after {attempts} copies")


def hadd_datasets_to_processes(
    era, input_dir, output_dir, add_derived_systs=True, dryRun=False,
    missing_only=False, selected_datasets=None,
):
    if not dryRun:
        import uproot
    selection = resolve_dataset_selection(ANALYSIS_PATH, era)

    if not dryRun and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Build the mapping Process -> list of associated datasets: the
    # process_names.yaml processes, plus the Combine processes that split one of
    # them (config/combine_processes.yaml, e.g. ST -> ST_s + ST_t).  Each file sums
    # the same keys (nominal and Up/Down) of its datasets, so merge-syst follows.
    process_mapping = dict(selection["process_datasets"])
    process_mapping.update(combine_datasets(process_mapping))
    # print(process_mapping)
    # 2. File check and execution/printing
    if dryRun:
        print("\n=== [DRY-RUN] MERGE PLAN (existing files only) ===")

    output_variants = (None, *dict.fromkeys(DY_COMPONENT_FILE_LABELS.values()))
    for process, datasets in process_mapping.items():
      for component_label in output_variants:
        datasets = [d for d in dict.fromkeys(datasets) if selected_datasets is None or d in selected_datasets]  # Remove duplicates

        # Do not duplicate the historical DY prefix in process-level files:
        # DY_DY_0J.root becomes DY_0J.root, while process-specific outputs use
        # e.g. DYto2Mu_MLL105To160_0J.root.
        output_component = (
            component_label.removeprefix("DY_") if component_label else None
        )
        output_suffix = f"_{output_component}" if output_component else ""
        output_name = f"{process}{output_suffix}.root"
        output_file_path = os.path.join(output_dir, output_name)
        if missing_only and os.path.exists(output_file_path):
            continue

        # Filter the datasets, keeping only those THAT EXIST ON DISK
        valid_dataset_files = []
        for dataset in datasets:
            for dataset_file_path in dataset_file_candidates(
                input_dir, process, dataset, component_label
            ):
                if os.path.exists(dataset_file_path):
                    valid_dataset_files.append(dataset_file_path)
                    break
        # IF THERE ARE NONE, PASS (skip the process entirely)
        if not valid_dataset_files:
            continue

        # --- IF DRY-RUN: only print what it would do ---
        if dryRun:
            print(f"\n📦 Expected output file: {output_name}")
            print(
                f"   ↳ Merge (hadd) of {len(valid_dataset_files)}/{len(datasets)} files found:"
            )
            for f in valid_dataset_files:
                print(f"     - {f}")
            continue

        # --- IF REAL EXECUTION ---
        print(
            f"📦 Creating {output_name} from {len(valid_dataset_files)} files..."
        )

        if len(valid_dataset_files) == 1:
            shutil.copyfile(valid_dataset_files[0], output_file_path)
            print(f"   -> Copied directly (single existing file).")
            continue

        try:
            histo_keys = []
            seen_keys = set()
            for file_path in valid_dataset_files:
                with uproot.open(file_path) as current_file:
                    for key in current_file.keys(recursive=True):
                        clean_key = key.split(";")[0]
                        if clean_key in seen_keys:
                            continue
                        try:
                            histo = current_file[clean_key]
                        except Exception:
                            continue
                        if hasattr(histo, "to_numpy") or "TH" in str(type(histo)):
                            histo_keys.append(clean_key)
                            seen_keys.add(clean_key)

            summed_histograms = {}

            for file_path in valid_dataset_files:
                with uproot.open(file_path) as current_file:
                    for key in histo_keys:
                        try:
                            histo = current_file[key]
                        except Exception:
                            continue
                        if hasattr(histo, "to_numpy") or "TH" in str(type(histo)):
                            current_hist = histo.to_hist()
                            if key not in summed_histograms:
                                summed_histograms[key] = current_hist
                            elif hist_sum_value(summed_histograms[key]) == 0.0 and hist_sum_value(current_hist) != 0.0:
                                summed_histograms[key] = current_hist
                            else:
                                summed_histograms[key] = summed_histograms[key] + current_hist

            write_histograms(output_file_path, summed_histograms)

            print(
                f"   ✅ Success! Saved to: {output_file_path}"
            )

        except Exception as e:
            print(
                f"   ❌ Error while processing process {process}: {e}"
            )

    if dryRun:
        print("\n===============================================================\n")
    else:
        if add_derived_systs:
            add_derived_systematics(era, output_dir)
        print("\n--- HADDing completed! ---")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( "--era", required=True, type=str)
    parser.add_argument( "--input-dir", required=True, type=str, help="ROOT file or dataset directory")
    parser.add_argument( "--output-dir", required=True, type=str, help="ROOT file or dataset directory")
    parser.add_argument( "--dryRun", action="store_true", help="dryRun only")
    parser.add_argument( "--add-derived-systs", action="store_true", help="add EWKZ unc")
    parser.add_argument(
        "--missing-only", action="store_true",
        help="create only process files that do not already exist",
    )
    parser.add_argument("--datasets", help="Comma-separated exact dataset whitelist for this campaign")
    args = parser.parse_args()

    hadd_datasets_to_processes(
        args.era, args.input_dir, args.output_dir, args.add_derived_systs,
        args.dryRun, args.missing_only,
        set(args.datasets.split(",")) if args.datasets else None,
    )
