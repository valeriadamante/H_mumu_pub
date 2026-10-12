#!/usr/bin/env python3
"""Distributions and ratios TAMU / Pisa for the VBF Z_CR_H_SB sync with Filippo, data or MC.

TAMU: our merged special sync tuple (VBF && Z_CR_H_SB); Pisa: Filippo's diagnostic
ntuple (first_failed_step == 0). The ratio uncertainty counts only the events selected
by one side. MC: --weight genWeight weights our events with genWeight and Filippo's
with --filippo-weight (a branch of his ntuple).

  python3 sync/python/plot_filippo_comparison.py --tag filippo_20260925 --output-dir <plots dir>
  python3 sync/python/plot_filippo_comparison.py --tag filippo_20260925 --dataset GluGluHto2Mu \
      --filippo-input '<his MC file pattern with {ds}>' --weight genWeight --filippo-weight genWeight \
      --output-dir <plots dir>
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_filippo_diagnostic import DATASETS, add_input_arguments, filippo_cut, filippo_steps, inputs, keys  # noqa: E402

PI = np.pi
# name, our branch, Filippo branch, bin edges, x label, overflow in last bin
VARIABLES = [
    ("m_mumu", "dimuon_mass", "Dimuon_mass", np.arange(76, 151, 2), r"$m_{\mu\mu}$ [GeV]", False),
    ("mu1_pt", "mu1_pt", "LeadingMuon_pt", np.arange(26, 151, 5), r"$p_T(\mu_1)$ [GeV]", True),
    ("mu2_pt", "mu2_pt", "SubLeadingMuon_pt", np.arange(20, 101, 4), r"$p_T(\mu_2)$ [GeV]", True),
    ("mu1_pt_corr", "mu1_pt_FSR_corr", "LeadingMuon_ptToUse", np.arange(20, 151, 5), r"corrected $p_T(\mu_1)$ [GeV]", True),
    ("mu2_pt_corr", "mu2_pt_FSR_corr", "SubLeadingMuon_ptToUse", np.arange(15, 101, 4), r"corrected $p_T(\mu_2)$ [GeV]", True),
    ("mu1_eta", "mu1_eta", "LeadingMuon_eta", np.linspace(-2.4, 2.4, 25), r"$\eta(\mu_1)$", False),
    ("mu2_eta", "mu2_eta", "SubLeadingMuon_eta", np.linspace(-2.4, 2.4, 25), r"$\eta(\mu_2)$", False),
    ("mu1_phi", "mu1_phi", "LeadingMuon_phi", np.linspace(-PI, PI, 25), r"$\phi(\mu_1)$", False),
    ("mu2_phi", "mu2_phi", "SubLeadingMuon_phi", np.linspace(-PI, PI, 25), r"$\phi(\mu_2)$", False),
    ("j1_pt", "leadingjet_pt", "LeadingJet_ptToUse", np.arange(35, 301, 10), r"$p_T(j_1)$ [GeV]", True),
    ("j2_pt", "subleadingjet_pt", "SubLeadingJet_ptToUse", np.arange(25, 151, 5), r"$p_T(j_2)$ [GeV]", True),
    ("j1_eta", "leadingjet_eta", "LeadingJet_eta", np.linspace(-4.7, 4.7, 48), r"$\eta(j_1)$", False),
    ("j2_eta", "subleadingjet_eta", "SubLeadingJet_eta", np.linspace(-4.7, 4.7, 48), r"$\eta(j_2)$", False),
    ("j1_phi", "leadingjet_phi", "LeadingJet_phi", np.linspace(-PI, PI, 25), r"$\phi(j_1)$", False),
    ("j2_phi", "subleadingjet_phi", "SubLeadingJet_phi", np.linspace(-PI, PI, 25), r"$\phi(j_2)$", False),
    ("mjj", "dijet_ls_mass", "Dijet_mass", np.arange(400, 2001, 50), r"$m_{jj}$ [GeV]", True),
    ("deta_jj", None, None, np.linspace(2.5, 8.5, 25), r"$|\Delta\eta_{jj}|$", False),
]
KEYS = ["run", "luminosityBlock", "event"]
TAMU, PISA = "#2a78d6", "#eb6834"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tag", default="filippo_20260925")
    parser.add_argument("--output-dir", type=Path, required=True)
    add_input_arguments(parser)
    parser.add_argument("--weight", choices=("none", "genWeight", "analysis"), default="none",
                        help="MC: weight our events with genWeight, or with the analysis weight of "
                             "--ours-weights (sync_mc_weights.py) (default: event counts)")
    parser.add_argument("--ours-weights", help="with --weight analysis: friend tuple of sync_mc_weights.py")
    parser.add_argument("--filippo-weight", help="MC: branch of Filippo's ntuple used as event weight with --weight")
    parser.add_argument("--title", help="legend title (default: the datasets and VBF && Z_CR_H_SB)")
    args = parser.parse_args()
    args.dataset = args.dataset or list(DATASETS)
    return args


def xrootd(path):
    """/eos/... -> root://: via FUSE uproot fa letture piccole e lente (28 s per
    un ramo delle nostre tuple contro 2 s via XRootD)."""
    path = str(path)
    for fuse, user in (("/eos/home-", "/eos/user/"),):
        if path.startswith(fuse):
            path = user + path[len(fuse):]
    if path.startswith("/eos/user/"):
        return "root://eosuser.cern.ch/" + path
    if path.startswith("/eos/cms/"):
        return "root://eoscms.cern.ch/" + path
    return path


def read(path, branches, cut):
    """Solo gli eventi che passano `cut`, a blocchi: i file di Pisa hanno 46 M
    eventi di cui ne passano 11 k, e leggere tutto in una volta satura memoria
    e XRootD."""
    url = xrootd(path)
    options = dict(handler=uproot.source.xrootd.XRootDSource, timeout=600) if url.startswith("root://") else {}
    with uproot.open(url, **options) as handle:
        parts = list(handle["Events"].iterate(branches, cut=cut, library="np", step_size="200 MB"))
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def load(args, cache):
    """Legge le tuple una volta sola e tiene gli eventi selezionati in `cache`."""
    if cache.is_file():
        with np.load(cache) as stored:
            arrays = {k: stored[k] for k in stored.files}
        tamu = {k[5:]: v for k, v in arrays.items() if k.startswith("tamu_")}
        pisa = {k[5:]: v for k, v in arrays.items() if k.startswith("pisa_")}
        print(f"inputs from the cache {cache} (delete it to re-read the tuples)")
    else:
        weighted = args.weight != "none"
        ours = [b for _, b, _, _, _, _ in VARIABLES if b] + (["genWeight"] if weighted else [])
        his = [b for _, _, b, _, _, _ in VARIABLES if b] + ["LeadingJet_eta", "SubLeadingJet_eta"] + \
            ([args.filippo_weight] if weighted and args.filippo_weight else [])
        tamu, pisa = [], []
        for ds in args.dataset:
            ours_path, fname = inputs(args, ds)
            tamu.append(read(ours_path, KEYS + ours, f"VBF & {args.region}"))
            with uproot.open(fname) as handle:
                names = filippo_steps(handle["Events"])
            pisa.append(read(fname, KEYS + his, filippo_cut(names, args.filippo_step)))
            print(f"{ds}: TAMU {len(tamu[-1]['run'])}, Pisa {len(pisa[-1]['run'])} selected events")
        tamu = {k: np.concatenate([a[k] for a in tamu]) for k in tamu[0]}
        pisa = {k: np.concatenate([a[k] for a in pisa]) for k in pisa[0]}
        np.savez(cache, **{f"tamu_{k}": v for k, v in tamu.items()},
                 **{f"pisa_{k}": v for k, v in pisa.items()})
    weighted = args.weight != "none"
    tamu["w"] = tamu["genWeight"] if weighted and "genWeight" in tamu else np.ones(len(tamu["run"]))
    if args.weight == "analysis":
        tamu["w"] = friend_weight(args.ours_weights, tamu)
    pisa["w"] = pisa[args.filippo_weight] if weighted and args.filippo_weight in pisa else np.ones(len(pisa["run"]))
    tamu["deta"] = np.abs(tamu["leadingjet_eta"] - tamu["subleadingjet_eta"])
    pisa["deta"] = np.abs(pisa["LeadingJet_eta"] - pisa["SubLeadingJet_eta"])
    kt = keys(*(tamu[k] for k in KEYS))
    kp = keys(*(pisa[k] for k in KEYS))
    tamu["only"] = ~np.isin(kt, kp)
    pisa["only"] = ~np.isin(kp, kt)
    return tamu, pisa


def friend_weight(friend_path, tamu, branch="weight_analysis"):
    """weight_analysis of the friend tuple, matched on run:lumi:event."""
    with uproot.open(friend_path) as handle:
        friend = handle["Events"].arrays(KEYS + [branch], library="np")
    fk = keys(*(friend[k] for k in KEYS))
    order = np.argsort(fk)
    tk = keys(*(tamu[k] for k in KEYS))
    pos = np.searchsorted(fk[order], tk)
    ok = (pos < len(fk)) & (fk[order][np.minimum(pos, len(fk) - 1)] == tk)
    if not ok.all():
        raise SystemExit(f"{(~ok).sum()} events without weight in {friend_path}")
    return friend[branch][order][pos]


def plot(tamu, pisa, name, ours, his, bins, label, overflow, output_dir, title):
    xt = np.clip(tamu[ours or "deta"], bins[0], bins[-1] - 1e-6)
    xp = np.clip(pisa[his or "deta"], bins[0], bins[-1] - 1e-6)
    wt, wp = tamu["w"], pisa["w"]
    ht, _ = np.histogram(xt, bins, weights=wt)
    hp, _ = np.histogram(xp, bins, weights=wp)
    stat_p, _ = np.histogram(xp, bins, weights=wp**2)
    only_t, _ = np.histogram(xt[tamu["only"]], bins, weights=wt[tamu["only"]]**2)
    only_p, _ = np.histogram(xp[pisa["only"]], bins, weights=wp[pisa["only"]]**2)
    centers, widths = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                      gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06))
    top.stairs(ht, bins, color=TAMU, lw=2, label=f"TAMU: {wt.sum():,.0f} events")
    top.errorbar(centers, hp, yerr=np.sqrt(stat_p), fmt="o", ms=4.5, color=PISA, mec="white", mew=0.6, lw=1.2,
                 label=f"Pisa: {wp.sum():,.0f} events")
    top.set_ylabel("Events / bin")
    top.set_ylim(0, max(ht.max(), hp.max()) * 1.45)
    top.grid(axis="y", color="#e6e6e3", lw=0.8)
    top.set_axisbelow(True)
    top.legend(frameon=False, loc="best", title=title, title_fontsize=10)

    ok = hp > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = ht / hp
        error = np.sqrt(only_t + only_p) / hp
    bottom.axhline(1, color="#8a8a85", lw=1)
    bottom.errorbar(centers[ok], ratio[ok], yerr=error[ok], xerr=widths[ok] / 2, fmt="o", ms=4, color="#1f1f1e", lw=1)
    deviation = np.max(np.abs(ratio[ok] - 1) + error[ok]) if ok.any() else 0.05
    span = max(0.05, np.ceil(deviation * 1.1 * 50) / 50)
    bottom.set_ylim(1 - span, 1 + span)
    bottom.set_ylabel("TAMU / Pisa")
    bottom.set_xlabel(label + ("  (last bin: overflow)" if overflow else ""))
    bottom.grid(axis="y", color="#e6e6e3", lw=0.8)
    bottom.text(0.01, 0.06, "uncertainty: non-shared events only", transform=bottom.transAxes, fontsize=8, color="#6b6b66")
    for ext in ("png", "pdf"):
        fig.savefig(output_dir / f"{name}.{ext}", bbox_inches="tight", dpi=130)
    plt.close(fig)


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
    cache = args.output_dir / f"inputs_{args.tag}_{'_'.join(args.dataset)}_{args.region}_{args.filippo_step}_{args.weight}.npz"
    tamu, pisa = load(args, cache)
    title = args.title or ("Run2024H" if sorted(args.dataset) == sorted(DATASETS) else ", ".join(args.dataset)) + \
        f", VBF && {args.region}"
    print(f"TAMU {len(tamu['run'])} events ({tamu['only'].sum()} only TAMU), "
          f"Pisa {len(pisa['run'])} events ({pisa['only'].sum()} only Pisa)")
    for variable in VARIABLES:
        plot(tamu, pisa, *variable, args.output_dir, title)
    print(f"wrote {len(VARIABLES)} plots (png + pdf) in {args.output_dir}")


if __name__ == "__main__":
    main()
