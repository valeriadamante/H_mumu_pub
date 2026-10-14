#!/usr/bin/env python3
"""Friend tuple with the analysis MC weights of an exchange tuple, for compare_filippo_weights.py.

The exchange tuples of the ``ggfvbf_pisa_20261013`` variant already carry the weights, so no
NanoAOD re-reading is needed (unlike sync_mc_weights.py).  The output has the same branch
names as the sync_mc_weights.py friend tuple.

  python3 sync/python/tuple_weights_friend.py TAMU-VBFHto2Mu_M125_powheg.root \\
      --selection "VBF & Signal_Fit" --output weights_VBF_Signal_Fit.root
"""
import argparse

import numpy as np
import uproot

KEYS = ["run", "luminosityBlock", "event"]
BRANCHES = {"genWeight": "weight_generator", "weight_pu": "weight_pileup",
            "weight_muIDIso": "weight_muon_idiso", "weight_trigger": "weight_trigger"}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tuple")
    parser.add_argument("--selection", default="VBF & Signal_Fit", help="uproot cut string")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    tree = uproot.open(args.tuple)["Events"]
    parts = list(tree.iterate(KEYS + list(BRANCHES.values()), cut=args.selection, library="np", step_size="200 MB"))
    data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    out = {k: data[k] for k in KEYS}
    out.update({friend: data[source].astype(np.float64) for friend, source in BRANCHES.items()})
    out["weight_analysis"] = out["genWeight"] * out["weight_pu"] * out["weight_muIDIso"] * out["weight_trigger"]
    with uproot.recreate(args.output) as handle:
        handle["Events"] = out
    print(f"{len(out['run'])} events, sum genWeight {out['genWeight'].sum():.6g}, sum weight_analysis {out['weight_analysis'].sum():.6g}")


if __name__ == "__main__":
    main()
