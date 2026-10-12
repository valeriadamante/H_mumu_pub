#!/usr/bin/env python3
"""Analysis MC weights for the events of a merged special sync tuple.

The special sync skim keeps only genWeight. This script adds, for the events
passing --selection, the weights of the Central entry of
config/Run3_2024/systematics.yaml, built with the framework's own functions:

  * weight_pu: corrections/pu.py table (floor(Pileup_nTrueInt), 2024 golden JSON).
    Pileup_nTrueInt is not in the tuple, so it is read from the NanoAOD files of
    the dataset (cached in --pileup-cache) and matched on run:lumi:event;
  * muon ID and isolation SFs: corrections/mu.py apply_muIDIso_weights, with
    the pT < 30 switch of systematics.yaml;
  * trigger SF: common/trigger_weights.py, single-leg IsoMu24 SF of the leading
    muon if it is trigger matched.

As in the production skim, mu{i}_pt and mu{i}_eta are the bsc + ScaRe + FSR
quantities (mu{i}_pt_FSR_corr, mu{i}_eta_FSR_corr); the SFs are evaluated at
mu{i}_pt_raw_noCorr. Writes a friend tuple (run, luminosityBlock, event and the
weights) and prints the yields.

  python3 sync/python/sync_mc_weights.py --dataset VBFHto2Mu_M125_powheg \\
      --output $HMUMU_SYNC/confronto_filippo_mc_vbfh/weights_VBF_Signal_Fit.root
"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))
os.environ.setdefault("ANALYSIS_PATH", str(REPOSITORY))
from common.paths import path  # noqa: E402

ERA = "Run3_2024"
KEYS = ["run", "luminosityBlock", "event"]
MUON_SF = ("((mu{i}_pt < 30) ? weight_mu{i}_TightPFIso_TightID_Central : weight_mu{i}_LoosePFIso_MediumID_Central)"
           " * ((mu{i}_pt < 30) ? weight_mu{i}_TightID_Trk_Central : weight_mu{i}_MediumID_Trk_Central)")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True, help="dataset of config/Run3_2024/samples_withfiles.yaml")
    parser.add_argument("--tag", default="filippo_20260925")
    parser.add_argument("--tuple", help="merged sync tuple (default: $HMUMU_SYNC/sync_tuple_<dataset>_<tag>.root)")
    parser.add_argument("--selection", default="VBF && Signal_Fit")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pileup-cache", type=Path,
                        help="npz with run, lumi, event, Pileup_nTrueInt of the NanoAOD (default: next to --output)")
    parser.add_argument("--workers", type=int, default=16)
    return parser.parse_args()


def read_pileup(url):
    import uproot
    for attempt in range(3):
        try:
            with uproot.open(url, handler=uproot.source.xrootd.XRootDSource, timeout=600) as handle:
                return handle["Events"].arrays(KEYS + ["Pileup_nTrueInt"], library="np")
        except Exception as error:  # transient XRootD errors
            if attempt == 2:
                raise RuntimeError(f"{url}: {error}") from error


def pileup_table(dataset, cache, workers):
    if cache.is_file():
        with np.load(cache) as stored:
            return {k: stored[k] for k in stored.files}
    import yaml
    files = yaml.safe_load((REPOSITORY / "config" / ERA / "samples_withfiles.yaml").read_text())[dataset]["filelist"]
    with ProcessPoolExecutor(workers) as pool:
        parts = list(pool.map(read_pileup, files))
    table = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    print(f"[weights] Pileup_nTrueInt of {len(table['run'])} NanoAOD events from {len(files)} files")
    np.savez(cache, **table)
    return table


def keys(run, lumi, event):
    run, lumi, event = (np.asarray(x).astype(np.uint64) for x in (run, lumi, event))
    return run << np.uint64(59) | lumi << np.uint64(39) | event


def main():
    args = arguments()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tuple_path = args.tuple or str(path("sync") / f"sync_tuple_{args.dataset}_{args.tag}.root")
    cache = args.pileup_cache or args.output.with_name(f"pileup_nanoaod_{args.dataset}.npz")
    pileup = pileup_table(args.dataset, cache, args.workers)

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gInterpreter.Declare(f'#include "{REPOSITORY / "analysis" / "AnalysisTools.h"}"')
    import common.utilities as utilities
    from corrections.mu import apply_muIDIso_weights
    from common.trigger_weights import AddTriggerWeightsAndErrors
    config = utilities.get_config(str(REPOSITORY / "config" / ERA / "maincfg.yaml"))

    df = ROOT.RDataFrame("Events", tuple_path).Filter(args.selection)
    for i in (1, 2):
        df = df.Redefine(f"mu{i}_pt", f"mu{i}_pt_FSR_corr").Redefine(f"mu{i}_eta", f"mu{i}_eta_FSR_corr")
    df, _ = apply_muIDIso_weights(df, config, False)
    df = AddTriggerWeightsAndErrors(df, False)
    df = df.Define("weight_muIDIso", " * ".join(MUON_SF.format(i=i) for i in (1, 2)))
    columns = KEYS + ["genWeight", "weight_muIDIso", "weight_TrgSF_singleMu_IsoMu24Central"]
    out = {k: np.asarray(v) for k, v in df.AsNumpy(columns).items()}
    out["weight_trigger"] = out.pop("weight_TrgSF_singleMu_IsoMu24Central")

    # PU: same table as corrections/pu.py, evaluated at floor(nTrueInt).
    import correctionlib
    from corrections.general import pog_folder_names, period_names
    period = period_names[config["era"]]
    pu_file = (f"/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/{pog_folder_names['LUM'][period]}"
               "/latest/puWeights_BCDEFGHI.json.gz")
    corr = correctionlib.CorrectionSet.from_file(pu_file)["Collisions24_BCDEFGHI_goldenJSON"]
    table = np.array([corr.evaluate(float(i), "nominal") for i in range(201)])
    nano_keys = keys(pileup["run"], pileup["luminosityBlock"], pileup["event"])
    order = np.argsort(nano_keys, kind="stable")
    ours = keys(out["run"], out["luminosityBlock"], out["event"])
    pos = np.searchsorted(nano_keys[order], ours)
    found = (pos < len(order)) & (nano_keys[order][np.minimum(pos, len(order) - 1)] == ours)
    if not found.all():
        raise SystemExit(f"{(~found).sum()} events of the tuple are not in the NanoAOD files")
    unique, counts = np.unique(nano_keys, return_counts=True)
    ambiguous = np.isin(ours, unique[counts > 1])
    n_true = pileup["Pileup_nTrueInt"][order][pos]
    out["Pileup_nTrueInt"] = n_true
    out["weight_pu"] = table[np.clip(np.floor(n_true).astype(int), 0, 200)]
    out["pu_ambiguous"] = ambiguous  # run:lumi:event reused by v2 and ext1: PU of one of the two
    out["weight_analysis"] = out["genWeight"] * out["weight_pu"] * out["weight_muIDIso"] * out["weight_trigger"]

    import uproot
    with uproot.recreate(args.output) as handle:
        handle["Events"] = out
    print(f"[weights] {len(ours)} events ({ambiguous.sum()} with a run:lumi:event reused in the NanoAOD), "
          f"wrote {args.output}")
    for name in ("genWeight", "weight_analysis"):
        print(f"  sum {name:16s} {out[name].sum():.6g}")
    for name in ("weight_pu", "weight_muIDIso", "weight_trigger"):
        print(f"  mean {name:15s} {out[name].mean():.5f}")


if __name__ == "__main__":
    main()
