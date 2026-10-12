#!/usr/bin/env python3
"""Per-event analysis weights of one Run3_2024 skim_v4 MC dataset, for the sync with gpizzati.

The weights are built by the framework (common.prepare_rdf.prepare_rdf), as in the histograms:
weight__Central = genWeight x lumi x xsec / sum of generator weights x PU x muon ID/iso SFs x
IsoMu24 trigger SF (x the DY reweights for DY, switched off with --no-dy-weights). The single
factors are saved too, so a disagreement can be traced to one of them:

  w_gen    genWeight
  w_pu     weight_pu
  w_idiso  muon ID x isolation SFs of both muons (pT < 30 GeV: tight, as in systematics.yaml)
  w_trg    weight_TrgSF_singleMu_IsoMu24Central
  w_norm   weight__Central / (w_gen w_pu w_idiso w_trg): lumi x xsec / sum(genWeight) (x DY weights)
  weight   weight__Central

Output: an .npz with run, luminosityBlock, event, the raw eta of the two muons (to match the
extensions) and the columns above, for every skim event (no selection).

  python3 sync/python/tamu_analysis_weights.py --dataset VBFHto2Mu_M125_powheg --output w.npz
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY))
os.environ.setdefault("ANALYSIS_PATH", str(REPOSITORY))
from common.paths import path  # noqa: E402

ERA = "Run3_2024"
IDISO = ("((mu{i}_pt < 30) ? weight_mu{i}_TightPFIso_TightID_Central : weight_mu{i}_LoosePFIso_MediumID_Central)"
         " * ((mu{i}_pt < 30) ? weight_mu{i}_TightID_Trk_Central : weight_mu{i}_MediumID_Trk_Central)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True, help="skim_v4 Run3_2024 dataset")
    parser.add_argument("--input-dir", type=Path, default=Path(str(path("skim"))) / ERA)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--file-list", type=Path, help="text file with the skim files to read (default: all); "
                        "the normalisation always uses the reports of the whole dataset")
    parser.add_argument("--no-dy-weights", action="store_true", help="DY: leave out the 012J, pT(ll) and N-jet weights")
    args = parser.parse_args()

    import ROOT
    from common.prepare_rdf import prepare_rdf
    from common.utilities import (get_config, get_segmentation_dict, initialize_root_runtime,
                                  list_root_files, report_path_for_root)
    initialize_root_runtime()
    ROOT.EnableImplicitMT(args.threads)
    config = REPOSITORY / "config" / ERA
    files = sorted(list_root_files(str(args.input_dir / args.dataset)))
    if not files:
        raise SystemExit(f"no ROOT file in {args.input_dir / args.dataset}")
    reports = [r for r in (report_path_for_root(f) for f in files) if r and os.path.isfile(r)]
    if len(reports) != len(files):
        raise SystemExit(f"{args.dataset}: {len(files) - len(reports)} files without report JSON")
    if args.file_list:
        files = args.file_list.read_text().split()
    dy = not args.no_dy_weights
    rdf = prepare_rdf(dataset_name=args.dataset, era=ERA,
                      selections_cfg=get_config(str(config / "selections.yaml")),
                      systematics_cfg=get_config(str(config / "systematics.yaml")),
                      is_data=False, input_files=files, seg_dict=get_segmentation_dict(reports),
                      skip_validation=True, enable_dy012j=dy, enable_dyptll=dy, enable_dynjets=dy)["inclusive"]
    rdf = (rdf.Define("w_gen", "double(genWeight)").Define("w_pu", "double(weight_pu)")
              .Define("w_idiso", f"double(({IDISO.format(i=1)}) * ({IDISO.format(i=2)}))")
              .Define("w_trg", "double(weight_TrgSF_singleMu_IsoMu24Central)")
              .Define("weight", "double(weight__Central)")
              .Define("w_norm", "w_gen * w_pu * w_idiso * w_trg != 0 ? weight / (w_gen * w_pu * w_idiso * w_trg) : 0."))
    columns = ["run", "luminosityBlock", "event", "mu1_eta_raw_noCorr", "mu2_eta_raw_noCorr",
               "w_gen", "w_pu", "w_idiso", "w_trg", "w_norm", "weight"]
    data = {k: np.asarray(v) for k, v in rdf.AsNumpy(columns).items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output, **data)
    print(f"{args.dataset}: {len(data['event'])} events, sum(weight) = {data['weight'].sum():.4g} -> {args.output}")


if __name__ == "__main__":
    main()
