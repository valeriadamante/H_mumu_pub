#!/usr/bin/env python3
"""Submit special_sync_skim.py to HTCondor, one job per NanoAOD file.

Outputs go to <output-root>/<dataset>/job_<n>.root and job_<n>_report.json.
Jobs whose report already exists on EOS are skipped, so rerunning the command
resubmits only the missing or failed ones.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import yaml

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))
from common.paths import path  # noqa: E402
DATASETS = ("Muon0_Run2024H", "Muon1_Run2024H")
GGFVBF_DATASETS = DATASETS + (
    "VBFHto2Mu_M125_powheg", "GluGluHto2Mu", "TTH_Hto2Mu",
    "ZH_Hto2Mu", "WplusH_Hto2Mu",
)


def requested_files(sample, exact_primary=False):
    """Return files for the requested DAS dataset, querying DAS if needed."""
    files = list(sample.get("filelist", []))
    datasets = sample.get("nanoAOD", [])
    if not exact_primary:
        return files
    if not datasets:
        raise RuntimeError("sample has no nanoAOD dataset")
    primary = datasets[0]
    # Existing lists are usable when they contain only the requested campaign.
    if files and len(datasets) == 1:
        return files
    if files:
        primary_processing = primary.split("/")[2]
        selected = [path for path in files if f"/{primary_processing}/" in path]
        if selected:
            return selected
    result = subprocess.run(
        ["dasgoclient", "-query", f"file dataset={primary}"],
        check=True, capture_output=True, text=True,
    )
    files = ["root://cms-xrd-global.cern.ch/" + line for line in result.stdout.splitlines() if line]
    if not files:
        raise RuntimeError(f"DAS returned no files for {primary}")
    return files


def xrootd(path):
    """Standard schedds refuse /eos paths in the arguments: use xrootd URLs."""
    path = path.replace("/eos/cms//store", "/eos/cms/store")
    if path.startswith("/eos/cms/"):
        return "root://eoscms.cern.ch/" + path
    if path.startswith("/eos/user/"):
        return "root://eosuser.cern.ch/" + path
    return path


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=("special_20260923", "filippo_20260925", "ggfvbf_20261015", "ggfvbf_pisa_20261013"),
                        default="special_20260923", help="selection of special_sync_skim.py")
    parser.add_argument("--tag", default=None, help="output folder name (default: the variant)")
    parser.add_argument("--output-root", type=Path, default=path("sync"))
    parser.add_argument("--group-name", default="vdamante",
                        help="prefix required by run3_sync: <GROUPNAME>-<DATASETNAME>")
    parser.add_argument("--proxy", default=str(REPOSITORY / "data" / "voms.proxy"))
    parser.add_argument("--flavour", default="workday")
    parser.add_argument("--max-jobs", type=int, default=None, help="submit at most N jobs (tests)")
    parser.add_argument("--no-submit", action="store_true")
    parser.add_argument("--dataset", action="append",
                        help="datasets of config/Run3_2024/samples_withfiles.yaml, data or MC; repeatable "
                             "(default: Muon0_Run2024H and Muon1_Run2024H)")
    parser.add_argument("--mc-lumi-split", action="store_true",
                        help="MC: apply the Run3_2024 event split of the production skim (see special_sync_skim.py)")
    return parser.parse_args()


def main():
    args = arguments()
    args.tag = args.tag or args.variant
    file_lists = REPOSITORY / "config" / "Run3_2024" / "samples_withfiles.yaml"
    if not file_lists.is_file():
        raise SystemExit(f"{file_lists} is missing: ./hmumu skim files --era Run3_2024 --use-ext")
    samples = yaml.safe_load(file_lists.read_text())
    log_dir = REPOSITORY / "htcondor" / "special_sync" / args.tag
    log_dir.mkdir(parents=True, exist_ok=True)
    lines = []
    unknown = [d for d in args.dataset or [] if d not in samples]
    if unknown:
        raise SystemExit(f"not in samples_withfiles.yaml: {unknown}")
    extra = "--mc-lumi-split" if args.mc_lumi_split else "none"
    default_datasets = GGFVBF_DATASETS if args.variant.startswith("ggfvbf") else DATASETS
    for dataset in args.dataset or default_datasets:
        output_dir = (args.output_root / f"{args.group_name}-{dataset}"
                      if args.variant.startswith("ggfvbf")
                      else args.output_root / args.tag / dataset)
        output_dir.mkdir(parents=True, exist_ok=True)
        files = requested_files(samples[dataset], exact_primary=args.variant.startswith("ggfvbf"))
        inputs_manifest = (output_dir / "inputs.txt" if args.variant.startswith("ggfvbf")
                           else output_dir.parent / f"{dataset}_inputs.txt")
        inputs_manifest.write_text(
            "".join(f"job_{i}\t{path}\n" for i, path in enumerate(files)))
        for index, path in enumerate(files):
            job = f"job_{index}"
            if (output_dir / f"{job}_report.json").is_file():
                continue
            lines.append(f"{args.proxy} {REPOSITORY} {dataset} {xrootd(path)} "
                         f"{xrootd(str(output_dir))} {job} {args.variant} {extra}")
    if args.max_jobs is not None:
        lines = lines[:args.max_jobs]
    # "input" is a Condor submit command (stdin transfer): never use it as a
    # queue variable, Condor would try to copy the NanoAOD file to the node.
    queue = log_dir / ("jobs_" + "_".join(args.dataset) + ".txt" if args.dataset else "jobs.txt")
    queue.write_text("\n".join(lines) + "\n")
    submit = queue.with_suffix(".sub")
    submit.write_text(f"""executable = {REPOSITORY}/sync/python/run_special_sync_skim.sh
arguments = $(proxy) $(repo) $(dataset) $(nanofile) $(outdir) $(job) $(variant) $(extra)
output = {log_dir}/$(dataset)_$(job).$(ClusterId).$(ProcId).out
error = {log_dir}/$(dataset)_$(job).$(ClusterId).$(ProcId).err
log = {log_dir}/special_sync.$(ClusterId).log
universe = vanilla
Requirements = (OpSysAndVer =?= "AlmaLinux9")
+JobFlavour = "{args.flavour}"
RequestCpus = 1
request_memory = 4GB
request_disk = 4GB
max_retries = 2
MY.SendCredential = true
batch_name = SpecialSync_{args.tag}
queue proxy, repo, dataset, nanofile, outdir, job, variant, extra from {queue}
""")
    print(f"[special-sync] {len(lines)} jobs, submit file {submit}")
    if lines and not args.no_submit:
        subprocess.run(["condor_submit", str(submit)], check=True)


if __name__ == "__main__":
    main()
