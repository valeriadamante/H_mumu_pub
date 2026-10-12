#!/usr/bin/env python3
"""Compare our analysis MC weights with Filippo's, factor by factor and in yield.

Ours: the friend tuple of sync_mc_weights.py (events of VBF && Signal_Fit).
Filippo: his diagnostic ntuple, events passing every step up to --filippo-step.
His per-step factors are matched to ours as
  BASESELECTION (genWeight)  <-> genWeight
  GoldenJson                 <-> weight_pu
  TwoMuons                   <-> weight_muIDIso (ID x isolation, both muons)
  LeadingMuonCut             <-> weight_trigger
On the events selected by both, the ratio ours / his of each factor is
summarised; the yields use each side's own events. Writes weights.md and
weight_ratios.csv in --output-dir.

  python3 sync/python/compare_filippo_weights.py \\
      --ours-weights $HMUMU_SYNC/confronto_filippo_mc_vbfh/weights_VBF_Signal_Fit.root \\
      --filippo-input '/eos/home-f/fcattafe/Hmm_sync_valeria/hmm_2024Summer24_VBF-Hto2Mu_M-125_powheg_*_VBF_SignalRegion_diagnostic_ntuples.root' \\
      --output-dir $HMUMU_SYNC/confronto_filippo_mc_vbfh
"""
import argparse
import csv
import glob
import sys
from pathlib import Path

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_filippo_diagnostic import filippo_cut, filippo_steps, keys  # noqa: E402

KEYS = ["run", "luminosityBlock", "event"]
PAIRS = [  # label, our branch, his branch
    ("genWeight", "genWeight", "BASESELECTION_weightFactor"),
    ("pileup", "weight_pu", "GoldenJson_weightFactor"),
    ("muon ID x iso (2 muons)", "weight_muIDIso", "TwoMuons_weightFactor"),
    ("trigger", "weight_trigger", "LeadingMuonCut_weightFactor"),
    ("total", "weight_analysis", None),
]


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ours-weights", required=True)
    parser.add_argument("--filippo-input", required=True, help="file or glob")
    parser.add_argument("--filippo-step", default="HiggsMass")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--xsec", type=float, default=0.00088819, help="pb, Filippo's (ours: 4.078 x 2.176e-4)")
    parser.add_argument("--lumi", type=float, default=109341.0, help="pb-1, Filippo's (ours: 109948.18)")
    return parser.parse_args()


def main():
    args = arguments()
    (filippo,) = glob.glob(args.filippo_input)
    with uproot.open(filippo) as handle:
        tree = handle["Events"]
        names = filippo_steps(tree)
        total_weight = f"{args.filippo_step}_cumulativeWeight"
        his = tree.arrays(KEYS + [b for _, _, b in PAIRS if b] + [total_weight],
                          cut=filippo_cut(names, args.filippo_step), library="np")
        sum_gen = float(handle["Runs"]["genEventSumw"].array(library="np").sum())
    his["total"] = his.pop(total_weight)
    with uproot.open(args.ours_weights) as handle:
        ours = handle["Events"].arrays(library="np")

    ko, kh = keys(*(ours[k] for k in KEYS)), keys(*(his[k] for k in KEYS))
    # run:lumi:event reused by v2 and ext1 cannot be matched one to one: drop them.
    def unique_mask(k):
        values, counts = np.unique(k, return_counts=True)
        return ~np.isin(k, values[counts > 1])
    uo, uh = unique_mask(ko), unique_mask(kh)
    common, io, ih = np.intersect1d(ko[uo], kh[uh], return_indices=True)
    io, ih = np.flatnonzero(uo)[io], np.flatnonzero(uh)[ih]

    rows = []
    for label, mine, theirs in PAIRS:
        a = ours[mine][io].astype(float)
        b = (his[theirs] if theirs else his["total"])[ih].astype(float)
        r = a / b
        rows.append((label, mine, theirs or total_weight, a.mean(), b.mean(), a.sum() / b.sum(),
                     *np.percentile(r, [1, 50, 99]), np.mean(np.abs(r - 1) < 1e-3)))

    norm = args.xsec * args.lumi / sum_gen
    yields = [
        ("TAMU, genWeight", ours["genWeight"].sum() * norm),
        ("TAMU, analysis weights", ours["weight_analysis"].sum() * norm),
        ("Pisa, genWeight", his["BASESELECTION_weightFactor"].sum() * norm),
        (f"Pisa, {total_weight}", his["total"].sum() * norm),
        ("TAMU, analysis weights, our xsec and lumi",
         ours["weight_analysis"].sum() * 4.078 * 2.176e-4 * 109948.18 / sum_gen),
    ]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "weight_ratios.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["factor", "ours", "his", "mean_ours", "mean_his", "sum_ratio",
                         "ratio_p01", "ratio_p50", "ratio_p99", "frac_within_0.1pc"])
        writer.writerows(rows)
    lines = [f"Events in common (unique run:lumi:event): {len(common)}; "
             f"sum genEventSumw = {sum_gen:.10g}; xsec x lumi = {args.xsec} pb x {args.lumi} pb-1.", "",
             "| factor | ours | Pisa | mean ours | mean Pisa | sum ours / sum Pisa | ratio p1 / p50 / p99 | within 0.1% |",
             "|---|---|---|---|---|---|---|---|"]
    for label, mine, theirs, ma, mb, sr, p1, p50, p99, frac in rows:
        lines.append(f"| {label} | `{mine}` | `{theirs}` | {ma:.4f} | {mb:.4f} | {sr:.4f} | "
                     f"{p1:.4f} / {p50:.4f} / {p99:.4f} | {100 * frac:.1f}% |")
    lines += ["", "| yield (VBF, 115-135 GeV) | events |", "|---|---|"]
    lines += [f"| {label} | {value:.3f} |" for label, value in yields]
    (args.output_dir / "weights.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
