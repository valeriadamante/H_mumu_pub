#!/usr/bin/env python3
"""Validate the files of a skim production.

The existence check used by the submitter is not enough: a job killed after the
``Snapshot`` leaves a full ``.root`` without ``report_*.json``, and a truncated
file opens but contains fewer events than those declared in the report.

For every chunk of the ``skim_chunks.json`` manifest it checks that

* the ``.root`` exists, is not empty and opens without errors;
* it contains the ``Events`` tree and the last event is readable;
* the ``report_*.json`` exists and is valid JSON;
* the tree entries match the last cut of the report.

Writes a JSON with the outcome of every chunk, reusable as a blacklist by
``tools/build_skim_overlay.py``.

    hmumu skim check --era 2024 \\
        --input /eos/.../skim_v4 --state-dir htcondor/skim_v4 \\
        --output results/skim_v4/validation_Run3_2024.json
"""

import argparse
import json
import os
from multiprocessing import Pool
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal


def check(task):
    """Return (dataset, index, outcome, entries, expected entries)."""
    directory, dataset, index = task
    root_path = os.path.join(directory, dataset, f"skim_{index}.root")
    report_path = os.path.join(directory, dataset, f"report_{index}.json")

    if not os.path.exists(root_path):
        return (dataset, index, "ROOT_MISSING", 0, -1)
    if os.path.getsize(root_path) == 0:
        return (dataset, index, "ROOT_EMPTY", 0, -1)

    try:
        tfile = ROOT.TFile.Open(root_path)
    except Exception as error:
        return (dataset, index, f"UNREADABLE ({type(error).__name__})", 0, -1)
    if not tfile or tfile.IsZombie():
        return (dataset, index, "ZOMBIE", 0, -1)

    tree = tfile.Get("Events")
    if not tree:
        tfile.Close()
        return (dataset, index, "NO_TREE_Events", 0, -1)

    entries = tree.GetEntries()
    try:
        tree.GetEntry(max(0, entries - 1))
    except Exception:
        tfile.Close()
        return (dataset, index, "LAST_EVENT_UNREADABLE", entries, -1)
    tfile.Close()

    if not os.path.exists(report_path):
        return (dataset, index, "JSON_MISSING", entries, -1)
    try:
        report = json.loads(Path(report_path).read_text())
    except Exception:
        return (dataset, index, "JSON_CORRUPT", entries, -1)

    passes = [v["pass"] for v in report.values() if isinstance(v, dict) and "pass" in v]
    expected = passes[-1] if passes else -1
    if expected >= 0 and expected != entries:
        return (dataset, index, f"MISMATCH report={expected} tree={entries}", entries, expected)

    return (dataset, index, "OK", entries, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--era", required=True)
    parser.add_argument("--input", required=True, type=Path, help="Base directory of the production.")
    parser.add_argument("--state-dir", required=True, type=Path,
                        help="State dir of the campaign, for the skim_chunks.json manifests.")
    parser.add_argument("--output", type=Path, help="JSON with the outcome of every chunk.")
    parser.add_argument("--datasets", nargs="*", help="Restrict to these datasets.")
    parser.add_argument("--jobs", type=int, default=16, help="Parallel processes.")
    parser.add_argument("--scan", action="store_true",
                        help="Ignore the manifests and validate all files present. Needed for "
                             "eras whose state dir covers only part of the datasets.")
    args = parser.parse_args()

    era_dir = args.input / args.era
    directory = str(era_dir)
    manifests = [] if args.scan else sorted(
        (args.state_dir / "log" / args.era).glob("*/skim_chunks.json"))
    tasks = []

    if manifests:
        for manifest in manifests:
            dataset = manifest.parent.name
            if args.datasets and dataset not in args.datasets:
                continue
            for chunk in json.loads(manifest.read_text())["chunks"]:
                tasks.append((directory, dataset, chunk["index"]))
    else:
        # No manifest: the era was produced with another state dir. Validate
        # what is there, without being able to notice chunks never written.
        if not args.scan:
            print(f"[WARN] no manifest under {args.state_dir}/log/{args.era}: "
                  f"validating the files present, chunks never produced cannot be detected")
        if not era_dir.is_dir():
            raise SystemExit(f"[ERROR] {era_dir} does not exist")
        for dataset_dir in sorted(p for p in era_dir.iterdir() if p.is_dir()):
            if args.datasets and dataset_dir.name not in args.datasets:
                continue
            for root_file in sorted(dataset_dir.glob("skim_*.root")):
                index = int(root_file.stem.split("_")[-1])
                tasks.append((directory, dataset_dir.name, index))

    if not tasks:
        raise SystemExit(f"[ERROR] no chunk to validate in {directory}")

    print(f"validating {len(tasks)} chunks of {args.era} in {directory}", flush=True)
    with Pool(args.jobs) as pool:
        results = pool.map(check, tasks, chunksize=4)

    bad = [r for r in results if r[2] != "OK"]
    entries = sum(r[3] for r in results if r[2] == "OK")

    print(f"\nchunks validated : {len(results)}")
    print(f"OK               : {len(results) - len(bad)}")
    print(f"PROBLEMS         : {len(bad)}")
    print(f"total events     : {entries}")

    if bad:
        print("\ndataset                                        idx  problem")
        for dataset, index, what, _, _ in bad:
            print(f"{dataset:45s} {index:4d}  {what}")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps([list(r) for r in results], indent=1))
        print(f"\noutcome written to {args.output}")

    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
