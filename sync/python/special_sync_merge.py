#!/usr/bin/env python3
"""Merge the special sync jobs into the final sync tuples and cutflows.

Refuses to merge unless every input file of both datasets has its report.
Writes, next to the other sync tuples:
  sync_tuple_<dataset>_<tag>.root (and sync_tuple_Run2024H_<tag>.root for the Muon0 + Muon1 data)
  <stem>_cutflow.csv, <stem>_cutflow_vbf_Z_CR_H_SB.csv, <stem>_categories_regions.csv
"""
import argparse
import csv
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))
from common.paths import path  # noqa: E402
from sync.python.prepare_tuple_forSync import write_cutflow  # noqa: E402

DATASETS = ("Muon0_Run2024H", "Muon1_Run2024H")
GGFVBF_DATASETS = DATASETS + (
    "VBFHto2Mu_M125_powheg", "GluGluHto2Mu", "TTH_Hto2Mu",
    "ZH_Hto2Mu", "WplusH_Hto2Mu",
)


def log(message):
    """Print a timestamped progress message immediately, also in batch logs."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] [special-sync] {message}", flush=True)


def human_size(size):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.1f} {unit}"
        size /= 1024


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="special_20260923")
    parser.add_argument("--output-root", type=Path, default=path("sync"))
    parser.add_argument("--group-name", default="vdamante")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dataset", action="append",
                        help="datasets to merge, data or MC; repeatable (default: Muon0_Run2024H and Muon1_Run2024H, "
                             "which also get the combined sync_tuple_Run2024H_<tag>)")
    return parser.parse_args()


def dataset_jobs(dataset_dir, inputs_path=None):
    manifest = inputs_path or dataset_dir / "inputs.txt"
    log(f"checking {dataset_dir.name} against {manifest}")
    inputs = manifest.read_text().split("\n")
    jobs = [line.split("\t")[0] for line in inputs if line]
    missing = [job for job in jobs if not (dataset_dir / f"{job}_report.json").is_file()
               or not (dataset_dir / f"{job}.root").is_file()]
    if missing:
        raise SystemExit(f"[ERROR] {dataset_dir.name}: {len(missing)}/{len(jobs)} jobs missing, e.g. {missing[:5]}")
    files = [dataset_dir / f"{job}.root" for job in jobs]
    log(f"{dataset_dir.name}: all {len(files)} job files and reports are present")
    return files


def summed(reports):
    report, cutflows, yields, total = {}, {}, {}, 0
    for path in reports:
        data = json.loads(path.read_text())
        for name, value in data["report"].items():
            report[name] = report.get(name, 0) + int(value if isinstance(value, int) else value["pass"])
        for key, steps in data["cutflows"].items():
            target = cutflows.setdefault(key, {})
            for step, count in steps:
                target[step] = target.get(step, 0) + count
        for kind, name, count in data["yields"]:
            yields[(kind, name)] = yields.get((kind, name), 0) + count
        total += data["events_after_skim"]
    return report, cutflows, yields, total


def entries(path):
    import ROOT
    handle = ROOT.TFile.Open(str(path))
    count = handle.Get("Events").GetEntries()
    handle.Close()
    return int(count)


def signature_groups(root_files):
    """Job files grouped by (branch, leaf type) signature; hadd silently drops files whose signature differs from the first one."""
    import ROOT
    groups = {}
    for path in root_files:
        handle = ROOT.TFile.Open(str(path))
        tree = handle.Get("Events")
        signature = frozenset((leaf.GetName(), leaf.GetTypeName()) for leaf in tree.GetListOfLeaves())
        handle.Close()
        groups.setdefault(signature, []).append(path)
    return groups


def unify_signatures(groups, output):
    """Snapshot every signature group onto the common columns, with sample_id as Long_t; returns the files to hadd."""
    import ROOT
    ROOT.EnableImplicitMT(4)
    common = set.intersection(*[{name for name, _ in signature} for signature in groups])
    first = next(iter(next(iter(groups.values()))))
    handle = ROOT.TFile.Open(str(first))
    columns = [leaf.GetName() for leaf in handle.Get("Events").GetListOfLeaves() if leaf.GetName() in common]
    handle.Close()
    dropped = set.union(*[{name for name, _ in signature} for signature in groups]) - common
    log(f"keeping {len(columns)} common columns, dropping {sorted(dropped)}")
    options = ROOT.RDF.RSnapshotOptions()
    options.fCompressionAlgorithm = ROOT.kLZMA
    options.fCompressionLevel = 1
    parts = []
    for index, files in enumerate(groups.values()):
        names = ROOT.std.vector("string")()
        for path in files:
            names.push_back(str(path))
        part = output.with_name(f"{output.stem}_unified_{index}.root")
        frame = ROOT.RDataFrame("Events", names)
        if "sample_id" in columns:
            frame = frame.Redefine("sample_id", "static_cast<Long_t>(sample_id)")
        frame.Snapshot("Events", str(part), columns, options)
        parts.append(part)
    return parts


def write_outputs(stem, root_files, overwrite):
    stem.parent.mkdir(parents=True, exist_ok=True)
    output = stem.with_suffix(".root")
    if output.exists() and not overwrite:
        raise SystemExit(f"[ERROR] {output} exists; pass --overwrite")
    original_files = list(root_files)
    input_size = sum(path.stat().st_size for path in root_files)
    log(f"merging {len(root_files)} ROOT files ({human_size(input_size)}) -> {output}")
    started = time.monotonic()
    groups = signature_groups(root_files)
    if len(groups) > 1:
        log(f"{len(groups)} different branch/type signatures among the job files: unifying them before hadd")
        root_files = unify_signatures(groups, output)
    process = subprocess.Popen(["hadd", "-f", str(output), *map(str, root_files)],
                               stdout=subprocess.DEVNULL)
    next_update = 30
    while process.poll() is None:
        time.sleep(1)
        elapsed = int(time.monotonic() - started)
        if elapsed >= next_update:
            partial_size = output.stat().st_size if output.exists() else 0
            log(f"hadd still running: {elapsed}s elapsed, output currently {human_size(partial_size)}")
            next_update += 30
    if process.returncode:
        raise subprocess.CalledProcessError(process.returncode, process.args)
    log(f"ROOT merge finished in {time.monotonic() - started:.1f}s ({human_size(output.stat().st_size)})")
    expected = sum(entries(path) for path in original_files)
    merged = entries(output)
    if merged != expected:
        raise SystemExit(f"[ERROR] {output}: {merged} entries, the job files hold {expected}; hadd dropped events")
    log(f"entries check passed: {merged} = sum over the job files")
    for path in root_files:
        if path not in original_files:
            path.unlink()
    reports = [path.with_name(f"{path.stem}_report.json") for path in original_files]
    log(f"summing {len(reports)} JSON reports and writing cutflows")
    report, cutflows, yields, total = summed(reports)
    report_counts = list(report.items())
    for key, steps in cutflows.items():
        suffix = "_cutflow" if key == "cutflow" else f"_{key}"
        write_cutflow(stem.with_name(stem.name + suffix + ".csv"), report_counts, list(steps.items()))
    with stem.with_name(stem.name + "_categories_regions.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["kind", "selection", "pass", "total_events", "efficiency"])
        for (kind, name), count in yields.items():
            writer.writerow([kind, name, count, total, count / total if total else 0.0])
    log(f"completed {output}")


def main():
    args = arguments()
    datasets = args.dataset or list(GGFVBF_DATASETS if args.tag.startswith("ggfvbf") else DATASETS)
    modern = args.tag.startswith("ggfvbf")
    log(f"starting tag={args.tag}, datasets={', '.join(datasets)}, overwrite={args.overwrite}")
    dirs = {dataset: (args.output_root / f"{args.group_name}-{dataset}" if modern
                      else args.output_root / args.tag / dataset) for dataset in datasets}
    per_dataset = {
        dataset: dataset_jobs(
            dirs[dataset], None if modern else args.output_root / args.tag / f"{dataset}_inputs.txt"
        ) for dataset in datasets
    }
    for dataset, files in per_dataset.items():
        log(f"processing dataset {dataset}")
        stem = dirs[dataset] / "sync_tuple" if modern else args.output_root / f"sync_tuple_{dataset}_{args.tag}"
        write_outputs(stem, files, args.overwrite)
    if set(DATASETS).issubset(datasets):
        log("processing combined Run2024H data tuple")
        stem = (args.output_root / f"{args.group_name}-Run2024H" / "sync_tuple" if modern
                else args.output_root / f"sync_tuple_Run2024H_{args.tag}")
        write_outputs(stem,
                      [f for dataset in DATASETS for f in per_dataset[dataset]], args.overwrite)
    vh = ("ZH_Hto2Mu", "WplusH_Hto2Mu")
    if set(vh).issubset(datasets):
        log("processing combined VH tuple")
        stem = (args.output_root / f"{args.group_name}-VH" / "sync_tuple" if modern
                else args.output_root / f"sync_tuple_VH_{args.tag}")
        write_outputs(stem,
                      [f for dataset in vh for f in per_dataset[dataset]], args.overwrite)
    log("all requested merges completed")


if __name__ == "__main__":
    main()
