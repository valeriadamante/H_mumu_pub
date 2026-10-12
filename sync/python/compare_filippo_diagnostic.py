#!/usr/bin/env python3
"""Compare a special sync tuple with Filippo's diagnostic ntuples.

Events are matched on run:luminosityBlock:event. Our selection is
`VBF && Z_CR_H_SB`; Filippo's is `first_failed_step == 0`. For the events
selected by only one side, the script prints where the other side loses them:
Filippo's first failed step, or our first failed step of cutflow_vbf_Z_CR_H_SB
(from the selection config stored in the job reports). For Filippo-only events
failing our electron veto it also counts the alternative electron vetoes stored
by the filippo_20260925 variant.

  python3 sync/python/compare_filippo_diagnostic.py --tag filippo_20260925 \
      --output-dir "$HMUMU_SYNC/filippo_comparison_20260925"
  # MC, once Filippo's ntuples exist ({ds} is replaced by each dataset):
  python3 sync/python/compare_filippo_diagnostic.py --tag filippo_20260925 --dataset GluGluHto2Mu \
      --filippo-input '/eos/home-f/fcattafe/Hmm_sync_valeria/<his MC file pattern with {ds}>' \
      --output-dir "$HMUMU_SYNC/filippo_comparison_mc"
  # MC signal region (Filippo's VBF_SignalRegion ntuple): our VBF && Signal_Fit
  # against his events passing every step up to HiggsMass (MC_selection is a
  # random ~50% event split that we do not reproduce):
  python3 sync/python/compare_filippo_diagnostic.py --dataset VBFHto2Mu_M125_powheg \
      --region Signal_Fit --filippo-step HiggsMass \
      --filippo-input '/eos/home-f/fcattafe/Hmm_sync_valeria/hmm_2024Summer24_VBF-Hto2Mu_M-125_powheg_*_VBF_SignalRegion_diagnostic_ntuples.root' \
      --output-dir "$HMUMU_SYNC/confronto_filippo_mc_vbfh"
"""
import argparse
import collections
import csv
import glob
import json
import sys
from pathlib import Path

import numpy as np
import uproot
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.paths import path  # noqa: E402

SYNC = path("sync")
FILIPPO = "/eos/home-f/fcattafe/Hmm_sync_valeria/hmm_2024Summer24_{ds}_*_diagnostic_ntuples.root"
DATASETS = ("Muon0_Run2024H", "Muon1_Run2024H")
OURS = "sync_tuple_{ds}_{tag}.root"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tag", default="filippo_20260925")
    parser.add_argument("--output-dir", type=Path, required=True)
    add_input_arguments(parser)
    return parser.parse_args()


def add_input_arguments(parser):
    """Inputs shared with plot_filippo_comparison.py; {ds} is replaced by each dataset."""
    parser.add_argument("--dataset", action="append",
                        help="dataset (special sync skim name), data or MC; repeatable (default: Muon0_Run2024H, Muon1_Run2024H)")
    parser.add_argument("--ours-tuple", default=str(SYNC / OURS),
                        help="our merged sync tuple; {ds} and {tag} are replaced (default: %(default)s)")
    parser.add_argument("--filippo-input", default=FILIPPO,
                        help="Filippo's diagnostic ntuple, file or glob; {ds} is replaced (default: %(default)s)")
    parser.add_argument("--region", default="Z_CR_H_SB",
                        help="our mass region, combined with VBF (default: %(default)s)")
    parser.add_argument("--filippo-step", default=None,
                        help="Filippo's events passing every step up to and including this one "
                             "(default: all steps, first_failed_step == 0)")


def inputs(args, ds):
    """(our tuple, Filippo's ntuple) for one dataset."""
    ours = Path(args.ours_tuple.format(ds=ds, tag=args.tag))
    found = sorted(glob.glob(args.filippo_input.format(ds=ds)))
    if len(found) != 1:
        raise SystemExit(f"{ds}: expected one Filippo ntuple for {args.filippo_input.format(ds=ds)}, found {found}")
    return ours, found[0]


RUN0 = 370000  # data, Run 3 2024: run - RUN0 < 2^15, lumi < 2^14, event < 2^34
MC_BIT = 63     # MC (run < 8): bit 63 set, run < 2^3, lumi < 2^20, event < 2^39


def keys(run, lumi, event):
    run, lumi, event = (x.astype(np.uint64) for x in (run, lumi, event))
    if len(run) and run.max() < 8:
        assert lumi.max() < 2**20 and event.max() < 2**39
        return np.uint64(1) << np.uint64(MC_BIT) | run << np.uint64(59) | lumi << np.uint64(39) | event
    assert run.min() >= RUN0 and (run - RUN0).max() < 2**15 and lumi.max() < 2**14 and event.max() < 2**34
    return (run - np.uint64(RUN0)) << np.uint64(48) | lumi << np.uint64(34) | event


def unpack(key):
    key = int(key)
    if key >> MC_BIT:
        return (key >> 59) & 7, (key >> 39) & (2**20 - 1), key & (2**39 - 1)
    return RUN0 + (key >> 48), (key >> 34) & (2**14 - 1), key & (2**34 - 1)


def filippo_steps(tree):
    """first_failed_step k is pass_* branch number k (0-based); BASESELECTION (0) always passes."""
    return [b.removeprefix("pass_") for b in tree.keys() if b.startswith("pass_")]


def filippo_pass(first_failed, names, step=None):
    """Events passing every step up to `step` (all steps if None)."""
    if step is None:
        return first_failed == 0
    return (first_failed == 0) | (first_failed > names.index(step))


def filippo_cut(names, step=None):
    """The same selection as an uproot cut string."""
    if step is None:
        return "first_failed_step == 0"
    return f"(first_failed_step == 0) | (first_failed_step > {names.index(step)})"


def our_steps(steps, region):
    """cutflow_vbf_Z_CR_H_SB with its last step (the mass region) replaced by `region`."""
    return [*steps[:-1], region] if steps[-1] != region else steps


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = args.dataset or list(DATASETS)
    report = json.loads(next((SYNC / args.tag / datasets[0]).glob("job_*_report.json")).read_text())
    steps = our_steps(yaml.safe_load(Path(report["selection_config"]).read_text())["cutflow_vbf_Z_CR_H_SB"],
                      args.region)
    electron_columns = ["n_goodElectrons", "n_goodElectrons_WP80", "n_goodElectrons_pt20", "n_goodElectrons_dRmu04"]
    summary, only_ours_rows, only_theirs_rows = [], [], []
    for ds in datasets:
        ours_path, filippo_path = inputs(args, ds)
        ours_tree = uproot.open(ours_path)["Events"]
        wanted = ["run", "luminosityBlock", "event", "VBF", args.region, "m_mumu", *steps, *electron_columns]
        ours = ours_tree.arrays([c for c in dict.fromkeys(wanted) if c in ours_tree.keys()], library="np")
        our_keys = keys(ours["run"], ours["luminosityBlock"], ours["event"])
        our_pass = ours["VBF"] & ours[args.region]

        tree = uproot.open(filippo_path)["Events"]
        names = filippo_steps(tree)
        theirs = tree.arrays(["run", "luminosityBlock", "event", "first_failed_step", "Dimuon_mass"], library="np")
        their_keys = keys(theirs["run"], theirs["luminosityBlock"], theirs["event"])
        their_pass = filippo_pass(theirs["first_failed_step"], names, args.filippo_step)
        for side, k, sel in (("ours", our_keys, our_pass), ("Filippo", their_keys, their_pass)):
            duplicates = int(sel.sum()) - len(np.unique(k[sel]))
            if duplicates:
                print(f"[WARNING] {ds}: {duplicates} duplicate run:lumi:event among the {side} selected events")

        a, b = set(our_keys[our_pass]), set(their_keys[their_pass])
        common, only_ours, only_theirs = a & b, a - b, b - a
        summary.append((ds, len(a), len(b), len(common), len(only_ours), len(only_theirs)))

        # Filippo's step for our extra events.
        order = np.argsort(their_keys)
        sorted_keys = their_keys[order]
        extra = np.array(sorted(only_ours), dtype=np.uint64)
        pos = np.searchsorted(sorted_keys, extra)
        found = (pos < len(sorted_keys)) & (sorted_keys[np.minimum(pos, len(sorted_keys) - 1)] == extra)
        failed = collections.Counter()
        for key, ok, p in zip(extra, found, pos):
            step = names[theirs["first_failed_step"][order[p]]] if ok else "non nella sua tupla"
            failed[step] += 1
            only_ours_rows.append((ds, *unpack(key), step))
        print(f"\n{ds}: solo nostri, suo primo passo fallito")
        for step, n in failed.most_common():
            print(f"  {step:30s} {n}")

        # Our first failed step for Filippo's extra events.
        index = {k: i for i, k in enumerate(our_keys)}
        lost = collections.Counter()
        electrons = collections.Counter()
        for key in only_theirs:
            i = index.get(key)
            if i is None:
                step = "non nella nostra tupla"
            else:
                step = next((s for s in steps if s in ours and not ours[s][i]), "?")
                if step == "electron_veto":
                    for column in electron_columns:
                        electrons[column] += int(ours[column][i] == 0)
            lost[step] += 1
            only_theirs_rows.append((ds, *unpack(key), step))
        print(f"{ds}: solo Filippo, nostro primo passo fallito")
        for step, n in lost.most_common():
            print(f"  {step:30s} {n}")
        if electrons:
            print(f"  veti elettroni alternativi che li salvano (su {lost['electron_veto']}): {dict(electrons)}")

    for name, rows, header in (("solo_nostri.csv", only_ours_rows, "suo_primo_passo_fallito"),
                               ("solo_filippo.csv", only_theirs_rows, "nostro_primo_passo_fallito")):
        with (args.output_dir / name).open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["dataset", "run", "luminosityBlock", "event", header])
            writer.writerows(rows)
    with (args.output_dir / "riassunto.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["dataset", "nostri", "Filippo", "comuni", "solo_nostri", "solo_Filippo"])
        writer.writerows(summary)
        total = ["+".join(datasets), *np.sum([row[1:] for row in summary], axis=0)]
        writer.writerow(total)
    print("\n" + "\n".join(",".join(map(str, row)) for row in [*summary, total]))


if __name__ == "__main__":
    main()
