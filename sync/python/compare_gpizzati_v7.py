#!/usr/bin/env python3
"""Sync of our skim_v4 (TAMU) with gpizzati hmm_skim v7, Run3_2024, data and MC.

The analysis selection of config/Run3_2024/selections.yaml is applied to both
skims, each with its own objects:

  step        TAMU (skim_v4 columns)                  gpizzati v7
  OS          muons_OS_and_presel                     mu1_charge * mu2_charge < 0
  trigger     muons_trigger (IsoMu24 + matching)      not stored: always true
  mu1_pt      mu1_pt > 26 (FSR_corr)                  mu1_pt > 26
  mu2_pt      mu2_pt > 20                             mu2_pt > 20
  ML          pfIsoId >= 2 && mediumId                same
  TT30        pT < 30 -> pfIsoId >= 4 && tightId      same
  bveto       SelectedJetTagSel                       < 1 medium and < 2 loose among Jet_good && !Jet_bad
  VBF_def     HasVBF, j1 >= 35, j2 >= 25              max-m_jj pair (>= 400, |deta| >= 2.5, 35/25) among Jet_good && !Jet_bad
  region      m_mumu                                  mll

Categories: baseline, ggF (= baseline && !VBF_def), VBF. Regions: Z_sideband,
H_sideband, Signal_Fit (data only with --unblind).

Only the lumisections present in both skims are compared, so a file missing on
one side does not look like a selection difference; the lumisections found in
one skim only are reported. Events are matched on run:luminosityBlock:event.

Outputs in <output-dir>/<dataset>/:
  README.md                 summary: coverage, yields, cutflows, agreement, where one-side events are lost
  yields.csv                per category and region: TAMU, gpizzati, common, only TAMU, only gpizzati
  cutflow.csv               cumulative counts step by step, both sides and common
  agreement.csv             on common selected events: fraction with the same value of each variable
  events_one_side_only.csv  every one-side event, the step where the other side loses it, and the
                            variables and step flags of both skims
  coverage.json             lumisections and events in both / one skim
  plots/<category>_<region>/*.png|pdf   only with --plots: distributions and ratio TAMU / gpizzati
                                        (uncertainty from the non-shared events only)

MC is unweighted unless --weight genWeight or --weight analysis. With --weight analysis both
sides carry all their analysis weights, normalised to the 2024 luminosity:
  TAMU      weight__Central from the framework (sync/python/tamu_analysis_weights.py, cached in
            <output-dir>/<dataset>/tamu_weights.npz): genWeight x lumi x xsec / sum(genWeight) x PU x
            muon ID/iso SFs x IsoMu24 trigger SF (x DY reweights for DY)
  gpizzati  genWeight x weight_sf_pu x weight_sf_mu{1,2}_{id,iso} x weight_trigger_SF x rescale
            (rescale = xsec x lumi / sumw from his results.json)
and the README compares the single factors (gen, pu, idiso, trg, norm) on the common events.

  source env.sh
  # one dataset, default inputs, summaries only:
  python3 sync/python/compare_gpizzati_v7.py --dataset Muon0_Run2024C \\
      --output-dir /eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7
  # explicit inputs (file, folder or glob) and plots:
  python3 sync/python/compare_gpizzati_v7.py --dataset VBFHto2Mu_M125_powheg \\
      --ours-input /eos/cms/store/group/phys_higgs/cmshmm/vdamante/skim_v4/Run3_2024/VBFHto2Mu_M125_powheg \\
      --theirs-input '/eos/cms/store/group/phys_higgs/cmshmm/gpizzati/hmm_skim/v7/2024/VBFHto2Mu_powheg_skimmed_*.root' \\
      --output-dir /eos/user/v/vdamante/H_mumu/Sync/confronto_gpizzati_v7 --plots
"""
import argparse
import glob
import json
import re
import subprocess
import sys
from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import uproot  # noqa: E402

OURS = "/eos/cms/store/group/phys_higgs/cmshmm/vdamante/skim_v4/Run3_2024/{ds}"
THEIRS = "/eos/cms/store/group/phys_higgs/cmshmm/gpizzati/hmm_skim/v7/2024"
KEYS = ["run", "luminosityBlock", "event"]
# MC extensions reuse run:lumi:event, so MC events are matched also on the raw eta of the two muons
MATCH = KEYS + ["eta_key"]
STEPS = ["OS", "trigger", "mu1_pt", "mu2_pt", "ML", "TT30", "bveto", "VBF_def", "region"]
REGIONS = {
    "Z_sideband": lambda m: (m > 70) & (m < 110),
    "H_sideband": lambda m: ((m > 110) & (m < 115)) | ((m > 135) & (m < 150)),
    "Signal_Fit": lambda m: (m > 115) & (m < 135),
}
CATEGORIES = ("baseline", "ggF", "VBF")
# our dataset name -> gpizzati v7 name, where they differ
NAMES = {
    "DYto2Mu_MLL_105to160_amcatnloFXFX": "DYTo2Mu_M_105to160_amcatnloFXFX",
    "DYto2Mu_MLL_105to160_amcatnloFXFX_Fil_VBF": "DYTo2Mu_M_105to160_amcatnloFXFX_vbf",
    "EWK_2Mu2J_MLL_105to160_herwig": "EWK_2L2J_M_105to160_madgraph_herwig",
    "GluGluHto2Mu": "GluGluHto2Mu_powheg",
    "VBFHto2Mu_M125_powheg": "VBFHto2Mu_powheg",
    "VBFHto2Mu_M125_amcatnlo": "VBFHto2Mu_amcatnlo",
    "WtoMuNu_amcatnloFXFX": "WToMuNu_amcatnloFXFX",
    "WtoENu_amcatnloFXFX": "WToENu_amcatnloFXFX",
    "WtoTauNu_amcatnloFXFX": "WToTauNu_amcatnloFXFX",
}
# variables compared on the common events: name, relative tolerance (0 = exact)
AGREEMENT = [("m_mumu", 1e-4), ("mu1_pt", 1e-4), ("mu2_pt", 1e-4), ("mu1_eta", 1e-4), ("n_jets", 0),
             ("j1_pt", 0.01), ("j2_pt", 0.01), ("mjj", 0.01)]
TAMU, GIORGIO = "#2a78d6", "#eb6834"
FACTORS = ("w_gen", "w_pu", "w_idiso", "w_trg", "w_norm")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", action="append", required=True,
                        help="our skim_v4 dataset name (Run3_2024), repeatable; e.g. Muon0_Run2024C, VBFHto2Mu_M125_powheg")
    parser.add_argument("--ours-input", help="TAMU input: file, folder or glob (default: skim_v4/Run3_2024/<dataset>); one --dataset only")
    parser.add_argument("--theirs-input", help="gpizzati input: file, folder or glob (default: v7/2024/<their name>_skimmed_*.root); one --dataset only")
    parser.add_argument("--output-dir", type=Path, required=True, help="summaries go to <output-dir>/<dataset>/")
    parser.add_argument("--plots", action="store_true", help="also draw the plots (default: CSV and markdown summaries only)")
    parser.add_argument("--plots-dir", type=Path, help="plots folder (default: <output-dir>/<dataset>/plots); implies --plots")
    kind = parser.add_mutually_exclusive_group()
    kind.add_argument("--data", action="store_true", help="treat the inputs as data (default: from the dataset name)")
    kind.add_argument("--mc", action="store_true", help="treat the inputs as MC")
    parser.add_argument("--max-files", type=int, help="read at most N files per side (tests; lumis are then restricted)")
    parser.add_argument("--weight", choices=("none", "genWeight", "analysis"), default="none", help="MC event weight for yields and plots")
    parser.add_argument("--no-dy-weights", action="store_true",
                        help="--weight analysis: TAMU weights without the DY reweights (gpizzati has none)")
    parser.add_argument("--threads", type=int, default=8, help="threads of the TAMU weight computation")
    parser.add_argument("--unblind", action="store_true", help="also compare data in Signal_Fit")
    parser.add_argument("--label-ours", default="TAMU (skim_v4)")
    parser.add_argument("--label-theirs", default="gpizzati (v7)")
    args = parser.parse_args()
    if (args.ours_input or args.theirs_input or args.plots_dir) and len(args.dataset) > 1:
        parser.error("--ours-input, --theirs-input and --plots-dir need a single --dataset")
    args.plots = args.plots or args.plots_dir is not None
    return args


def theirs_name(ds):
    """Data: Muon0_Run2024C -> Muon0_Run2024C_v1; Run2024I_v1/v2 keep their suffix."""
    if re.fullmatch(r"Muon[01]_Run2024[A-Z]", ds):
        return ds + "_v1"
    return NAMES.get(ds, ds)


def files(path, max_files, prefix=None):
    """ROOT files from a file, a folder (only <prefix>*.root if given and present) or a glob."""
    if Path(path).is_dir():
        found = glob.glob(f"{path}/{prefix}*.root") if prefix else []
        found = found or glob.glob(f"{path}/*.root")
    else:
        found = glob.glob(path)
    number = lambda f: int((re.findall(r"(\d+)\.root$", f) or ["0"])[0])  # noqa: E731
    found = sorted(found, key=lambda f: (number(f), f))
    if not found:
        raise FileNotFoundError(f"no ROOT file in {path}")
    return found[:max_files] if max_files else found


def input_totals(paths, side):
    """Events before and after the skim, from the JSON written next to each ROOT file."""
    count = lambda v: v if isinstance(v, (int, float)) else v.get("pass", 0)  # noqa: E731
    before = after_split = after = 0
    found = 0
    for path in paths:
        if side == "TAMU":
            report = Path(path).with_name(Path(path).name.replace("skim_", "report_")).with_suffix(".json")
            if report.exists():
                steps = json.loads(report.read_text())
                before += count(steps["Initial"])
                after_split += count(next((v for k, v in steps.items() if k.startswith("Orthogonal_lumi_split")), steps["Initial"]))
                after += count(list(steps.values())[-1])
                found += 1
        else:
            report = Path(path.replace(".root", "_output.json"))
            if report.exists():
                info = next(iter(json.loads(report.read_text()).values()))
                before += info["nevents_before"]
                after_split += info["nevents_before"]
                after += info["nevents_after"]
                found += 1
    return {"reports": found, "before_skim": before, "after_lumi_split": after_split, "after_skim": after}


def read(paths, branches):
    """Concatenate the Events trees; files without events (no branches) are skipped.
    EOS I/O errors are retried a few times before failing."""
    import time
    parts = []
    for path in paths:
        for attempt in range(4):
            try:
                tree = uproot.open(path)["Events"]
                if tree.num_entries and tree.keys():
                    parts.append(tree.arrays([b for b in branches if b in tree.keys()]))
                break
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(30 * (attempt + 1))
    return ak.concatenate(parts)


def pair_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
    def p4(pt, eta, phi, m):
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        return np.sqrt(px**2 + py**2 + pz**2 + m**2), px, py, pz
    e1, x1, y1, z1 = p4(pt1, eta1, phi1, m1)
    e2, x2, y2, z2 = p4(pt2, eta2, phi2, m2)
    return np.sqrt(np.maximum((e1 + e2)**2 - (x1 + x2)**2 - (y1 + y2)**2 - (z1 + z2)**2, 0))


def eta_key(eta1, eta2):
    """Order-independent integer made of the two raw muon eta (NanoAOD values, same on both sides)."""
    a, b = (np.round(ak.to_numpy(x) * 1000).astype(np.int64) for x in (eta1, eta2))
    return np.minimum(a, b) * 10000 + np.maximum(a, b)


MUON = [f"mu{i}_{v}" for i in (1, 2) for v in ("pt", "eta", "phi", "charge", "pfIsoId", "mediumId", "tightId")]


def common_columns(a, pt, eta, weight):
    """Muon and jet variables used by plots and diagnostics, identical names on both sides."""
    first = lambda x, i: ak.to_numpy(ak.fill_none(ak.pad_none(x, i + 1)[:, i], np.nan))  # noqa: E731
    out = {k: ak.to_numpy(a[k]) for k in KEYS}
    for i in (1, 2):
        for v in ("pt", "eta", "phi"):
            out[f"mu{i}_{v}"] = ak.to_numpy(a[f"mu{i}_{v}"])
    for i in (0, 1):
        out[f"j{i + 1}_pt"], out[f"j{i + 1}_eta"] = first(pt, i), first(eta, i)
    out["n_jets"] = ak.to_numpy(ak.num(pt))
    out["weight"] = weight
    return out


def load_ours(paths, weighted):
    jets = [f"SelectedJet_{v}" for v in ("pt", "eta", "phi", "mass")]
    a = read(paths, KEYS + MUON + jets + ["muons_OS_and_presel", "muons_trigger", "SelectedJetTagSel", "HasVBF",
                                          "VBFJetIdx_1", "VBFJetIdx_2", "m_mumu", "genWeight",
                                          "mu1_eta_raw_noCorr", "mu2_eta_raw_noCorr"])
    pt, eta, phi, mass = (a[j] for j in jets)
    i1, i2 = ak.to_numpy(a.VBFJetIdx_1), ak.to_numpy(a.VBFJetIdx_2)
    has = ak.to_numpy(a.HasVBF).astype(bool)
    padded = {v: ak.to_numpy(ak.fill_none(ak.pad_none(x, 20, clip=True), np.nan)) for v, x in
              (("pt", pt), ("eta", eta), ("phi", phi), ("mass", mass))}
    rows = np.arange(len(a))
    j1 = {v: np.where(has, padded[v][rows, np.clip(i1, 0, 19)], np.nan) for v in padded}
    j2 = {v: np.where(has, padded[v][rows, np.clip(i2, 0, 19)], np.nan) for v in padded}
    weight = ak.to_numpy(a.genWeight) if weighted and "genWeight" in a.fields else np.ones(len(a))
    out = common_columns(a, pt, eta, weight)
    out["m_mumu"] = ak.to_numpy(a.m_mumu)
    out["eta_key"] = eta_key(a.mu1_eta_raw_noCorr, a.mu2_eta_raw_noCorr)
    out["mjj"] = pair_mass(j1["pt"], j1["eta"], j1["phi"], j1["mass"], j2["pt"], j2["eta"], j2["phi"], j2["mass"])
    out["deta_jj"] = np.abs(j1["eta"] - j2["eta"])
    steps = {"OS": ak.to_numpy(a.muons_OS_and_presel), "trigger": ak.to_numpy(a.muons_trigger),
             "bveto": ak.to_numpy(a.SelectedJetTagSel), "VBF_def": has & (j1["pt"] >= 35) & (j2["pt"] >= 25)}
    return muon_steps(a, out, steps)


THEIR_WEIGHTS = ["weight_sf_pu", "weight_trigger_SF"] + [f"weight_sf_mu{i}_{k}" for i in (1, 2) for k in ("id", "iso")]


def load_theirs(paths, weighted):
    jets = [f"Jet_{v}" for v in ("pt", "eta", "phi", "mass", "good", "bad", "btag_L", "btag_M")]
    a = read(paths, KEYS + MUON + jets + ["mll", "genWeight", "mu1_eta_no_fsr", "mu2_eta_no_fsr"]
             + (THEIR_WEIGHTS if weighted == "analysis" else []))
    sel = a.Jet_good & ~a.Jet_bad  # already cleaned from the muons, pT ordered
    pt, eta, phi, mass = a.Jet_pt[sel], a.Jet_eta[sel], a.Jet_phi[sel], a.Jet_mass[sel]
    # production VBF pair: highest m_jj among the pairs passing the thresholds
    first, second = ak.unzip(ak.combinations(ak.local_index(pt), 2))
    p1 = (pt[first], eta[first], phi[first], mass[first])
    p2 = (pt[second], eta[second], phi[second], mass[second])
    mjj = pair_mass(*p1, *p2)
    ok = (mjj >= 400) & (abs(p1[1] - p2[1]) >= 2.5) & (p1[0] >= 35) & (p2[0] >= 25)
    best = ak.argmax(ak.where(ok, mjj, -1), axis=1, keepdims=True)
    has = ak.to_numpy(ak.fill_none(ak.firsts(ok[best]), False))
    pick = lambda x: ak.to_numpy(ak.fill_none(ak.firsts(x[best]), np.nan))  # noqa: E731
    weight = ak.to_numpy(a.genWeight) if weighted and "genWeight" in a.fields else np.ones(len(a))
    out = common_columns(a, pt, eta, weight)
    if weighted == "analysis":
        col = lambda name: ak.to_numpy(a[name]).astype(np.float64)  # noqa: E731
        out.update(w_gen=col("genWeight"), w_pu=col("weight_sf_pu"), w_trg=col("weight_trigger_SF"),
                   w_idiso=col("weight_sf_mu1_id") * col("weight_sf_mu1_iso") * col("weight_sf_mu2_id") * col("weight_sf_mu2_iso"))
    out["m_mumu"] = ak.to_numpy(a.mll)
    out["eta_key"] = eta_key(a.mu1_eta_no_fsr, a.mu2_eta_no_fsr)
    out["mjj"] = np.where(has, pick(mjj), np.nan)
    out["deta_jj"] = np.where(has, np.abs(pick(p1[1]) - pick(p2[1])), np.nan)
    central = abs(eta) < 2.5
    n_loose = ak.to_numpy(ak.sum(a.Jet_btag_L[sel] & central, axis=1))
    n_medium = ak.to_numpy(ak.sum(a.Jet_btag_M[sel] & central, axis=1))
    steps = {"OS": ak.to_numpy(a.mu1_charge * a.mu2_charge < 0), "trigger": np.ones(len(a), bool),
             "bveto": (n_medium < 1) & (n_loose < 2), "VBF_def": has}
    return muon_steps(a, out, steps)


def attach_tamu_weights(ours, ds, args, out_dir, ours_files):
    """Merge the framework analysis weights (tamu_analysis_weights.py, cached) on run:lumi:event:eta."""
    cache = out_dir / ("tamu_weights_nodyw.npz" if args.no_dy_weights else "tamu_weights.npz")
    if not cache.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(Path(__file__).with_name("tamu_analysis_weights.py")), "--dataset", ds,
                   "--output", str(cache), "--threads", str(args.threads)] + (["--no-dy-weights"] if args.no_dy_weights else [])
        if args.ours_input and Path(args.ours_input).is_dir():
            command += ["--input-dir", str(Path(args.ours_input).parent)]
        if args.max_files:
            listing = out_dir / "tamu_files.txt"
            listing.write_text("\n".join(ours_files) + "\n")
            command += ["--file-list", str(listing)]
        subprocess.run(command, check=True)
    w = np.load(cache)
    table = pd.DataFrame({k: w[k] for k in KEYS + list(FACTORS) + ["weight"]})
    table["eta_key"] = eta_key(w["mu1_eta_raw_noCorr"], w["mu2_eta_raw_noCorr"])
    table = table.drop_duplicates(MATCH, keep=False).set_index(MATCH)
    found = table.reindex(pd.MultiIndex.from_frame(ours[MATCH]))
    missing = int(found["weight"].isna().sum())
    for column in list(FACTORS) + ["weight"]:
        ours[column] = found[column].values
    return ours, missing


def theirs_rescale(ds, theirs_files):
    info = json.loads((Path(theirs_files[0]).parent / "results.json").read_text())
    return float(info[theirs_name(ds)]["rescale"])


def muon_steps(a, out, steps):
    pt1, pt2 = out["mu1_pt"], out["mu2_pt"]
    col = lambda name: ak.to_numpy(a[name])  # noqa: E731
    ml = lambda i: (col(f"mu{i}_pfIsoId") >= 2) & col(f"mu{i}_mediumId").astype(bool)  # noqa: E731
    tt = lambda i, pt: (pt >= 30) | ((col(f"mu{i}_pfIsoId") >= 4) & col(f"mu{i}_tightId").astype(bool))  # noqa: E731
    steps.update({"mu1_pt": pt1 > 26, "mu2_pt": pt2 > 20, "ML": ml(1) & ml(2), "TT30": tt(1, pt1) & tt(2, pt2)})
    table = pd.DataFrame(out)
    for step, value in steps.items():
        table[f"pass_{step}"] = np.asarray(value, bool)
    return table


def passes(table, category, region):
    """Ordered step -> boolean array for one category/region."""
    result = {s: table[f"pass_{s}"].values for s in STEPS if s not in ("VBF_def", "region")}
    if category == "VBF":
        result["VBF_def"] = table.pass_VBF_def.values
    elif category == "ggF":
        result["!VBF_def"] = ~table.pass_VBF_def.values
    result["region"] = REGIONS[region](table.m_mumu.values)
    return result


def first_failed(steps):
    result = np.full(len(next(iter(steps.values()))), "PASS", dtype=object)
    for step, ok in reversed(list(steps.items())):
        result[~ok] = step
    return result


VARIABLES = [
    ("m_mumu", {"Z_sideband": np.arange(70, 110.1, 1), "H_sideband": np.arange(110, 150.1, 2),
                "Signal_Fit": np.arange(115, 135.1, 1)}, r"$m_{\mu\mu}$ [GeV]", False),
    ("mu1_pt", np.arange(26, 151, 5), r"$p_T(\mu_1)$ [GeV]", True),
    ("mu2_pt", np.arange(20, 101, 4), r"$p_T(\mu_2)$ [GeV]", True),
    ("mu1_eta", np.linspace(-2.4, 2.4, 25), r"$\eta(\mu_1)$", False),
    ("mu2_eta", np.linspace(-2.4, 2.4, 25), r"$\eta(\mu_2)$", False),
    ("mu1_phi", np.linspace(-np.pi, np.pi, 25), r"$\phi(\mu_1)$", False),
    ("mu2_phi", np.linspace(-np.pi, np.pi, 25), r"$\phi(\mu_2)$", False),
    ("n_jets", np.arange(-0.5, 7.5, 1), r"$N_{jets}$", True),
    ("j1_pt", np.arange(25, 301, 10), r"$p_T(j_1)$ [GeV]", True),
    ("j2_pt", np.arange(25, 151, 5), r"$p_T(j_2)$ [GeV]", True),
    ("j1_eta", np.linspace(-4.7, 4.7, 48), r"$\eta(j_1)$", False),
    ("j2_eta", np.linspace(-4.7, 4.7, 48), r"$\eta(j_2)$", False),
    ("mjj", np.arange(400, 2001, 50), r"$m_{jj}$ of the VBF pair [GeV]", True),
    ("deta_jj", np.linspace(2.5, 8.5, 25), r"$|\Delta\eta_{jj}|$ of the VBF pair", False),
]


def plot(ours, theirs, only_ours, only_theirs, category, region, labels, directory):
    directory.mkdir(parents=True, exist_ok=True)
    for name, bins, xlabel, overflow in VARIABLES:
        if isinstance(bins, dict):
            bins = bins[region]
        if name in ("mjj", "deta_jj") and category != "VBF":
            continue
        xo, xt = ours[name].values, theirs[name].values
        wo, wt = ours.weight.values, theirs.weight.values
        fo, ft = np.isfinite(xo), np.isfinite(xt)
        if not fo.any() and not ft.any():
            continue
        clip = lambda x: np.clip(x, bins[0], bins[-1] - 1e-6)  # noqa: E731
        ho, _ = np.histogram(clip(xo[fo]), bins, weights=wo[fo])
        ht, _ = np.histogram(clip(xt[ft]), bins, weights=wt[ft])
        so, _ = np.histogram(clip(xo[fo & only_ours]), bins, weights=wo[fo & only_ours]**2)
        st, _ = np.histogram(clip(xt[ft & only_theirs]), bins, weights=wt[ft & only_theirs]**2)
        stat_t, _ = np.histogram(clip(xt[ft]), bins, weights=wt[ft]**2)
        centers, widths = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)

        fig, (top, bottom) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                          gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06))
        top.stairs(ho, bins, color=TAMU, lw=2, label=f"{labels[0]}: {wo[fo].sum():,.0f}")
        top.errorbar(centers, ht, yerr=np.sqrt(stat_t), fmt="o", ms=4.5, color=GIORGIO, mec="white", mew=0.6,
                     lw=1.2, label=f"{labels[1]}: {wt[ft].sum():,.0f}")
        top.set_ylabel("Events / bin")
        top.set_ylim(0, max(ho.max(), ht.max(), 1) * 1.45)
        top.grid(axis="y", color="#e6e6e3", lw=0.8)
        top.set_axisbelow(True)
        top.legend(frameon=False, loc="best", title=f"Run3_2024, {category}, {region}", title_fontsize=10)
        ok = ht > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio, error = ho / ht, np.sqrt(so + st) / ht
        bottom.axhline(1, color="#8a8a85", lw=1)
        bottom.errorbar(centers[ok], ratio[ok], yerr=error[ok], xerr=widths[ok] / 2, fmt="o", ms=4, color="#1f1f1e", lw=1)
        deviation = np.nanmax(np.abs(ratio[ok] - 1) + error[ok]) if ok.any() else 0.05
        span = min(max(0.05, np.ceil(deviation * 1.1 * 50) / 50), 1.0)
        bottom.set_ylim(1 - span, 1 + span)
        bottom.set_ylabel("TAMU / gpizzati")
        bottom.set_xlabel(xlabel + ("  (last bin: overflow)" if overflow else ""))
        bottom.grid(axis="y", color="#e6e6e3", lw=0.8)
        bottom.text(0.01, 0.06, "uncertainty: non-shared events only", transform=bottom.transAxes, fontsize=8, color="#6b6b66")
        for ext in ("png", "pdf"):
            fig.savefig(directory / f"{name}.{ext}", bbox_inches="tight", dpi=130)
        plt.close(fig)


def compare(ds, args):
    is_data = args.data or (not args.mc and ds.startswith("Muon"))
    weighted = (args.weight if args.weight != "none" else False) if not is_data else False
    ours_files = files(args.ours_input or OURS.format(ds=ds), args.max_files)
    theirs_files = files(args.theirs_input or THEIRS, args.max_files, prefix=f"{theirs_name(ds)}_skimmed_")
    totals = {"TAMU": input_totals(ours_files, "TAMU"), "gpizzati": input_totals(theirs_files, "gpizzati")}
    ours, theirs = load_ours(ours_files, weighted), load_theirs(theirs_files, weighted)
    missing_weights = 0
    if weighted == "analysis":
        ours, missing_weights = attach_tamu_weights(ours, ds, args, args.output_dir / ds, ours_files)
        theirs["w_norm"] = theirs_rescale(ds, theirs_files)
        theirs["weight"] = theirs[list(FACTORS)].prod(axis=1)
        ours[["weight", *FACTORS]] = ours[["weight", *FACTORS]].fillna(0.0)
    if is_data:
        ours["eta_key"], theirs["eta_key"] = 0, 0
    # MC samples with extensions reuse run:lumi:event for different events: those keys cannot be
    # matched, so they are dropped on both sides and counted. In data a duplicate is an error.
    duplicated = {}
    for name, t in (("TAMU", ours), ("gpizzati", theirs)):
        mask = t.duplicated(MATCH, keep=False)
        if mask.any() and is_data:
            raise RuntimeError(f"{ds}: duplicated events in the {name} skim")
        duplicated[name] = pd.MultiIndex.from_frame(t.loc[mask, MATCH])
    bad = duplicated["TAMU"].union(duplicated["gpizzati"])
    dropped = {}
    for name, t in (("TAMU", ours), ("gpizzati", theirs)):
        mask = pd.MultiIndex.from_frame(t[MATCH]).isin(bad)
        dropped[name] = int(mask.sum())
        if name == "TAMU":
            ours = t[~mask].reset_index(drop=True)
        else:
            theirs = t[~mask].reset_index(drop=True)
    lumis_o = set(zip(ours.run, ours.luminosityBlock))
    lumis_t = set(zip(theirs.run, theirs.luminosityBlock))
    shared = lumis_o & lumis_t
    in_shared = lambda t: np.fromiter((k in shared for k in zip(t.run, t.luminosityBlock)), bool, len(t))  # noqa: E731
    keep_o, keep_t = in_shared(ours), in_shared(theirs)
    coverage = {"is_data": bool(is_data), "weights": weighted or "none",
                "TAMU_events_without_analysis_weight": missing_weights,
                "TAMU_input": args.ours_input or OURS.format(ds=ds), "TAMU_files": len(ours_files),
                "gpizzati_input": args.theirs_input or f"{THEIRS}/{theirs_name(ds)}_skimmed_*.root",
                "gpizzati_files": len(theirs_files), "input_totals": totals,
                "lumis_both": len(shared), "lumis_only_TAMU": len(lumis_o - shared), "lumis_only_gpizzati": len(lumis_t - shared),
                "events_TAMU_outside_shared_lumis": int((~keep_o).sum()),
                "events_gpizzati_outside_shared_lumis": int((~keep_t).sum()),
                "events_TAMU_dropped_duplicate_keys": dropped["TAMU"],
                "events_gpizzati_dropped_duplicate_keys": dropped["gpizzati"]}
    ours, theirs = ours[keep_o].reset_index(drop=True), theirs[keep_t].reset_index(drop=True)
    index_o, index_t = pd.MultiIndex.from_frame(ours[MATCH]), pd.MultiIndex.from_frame(theirs[MATCH])
    in_t, in_o = index_o.isin(index_t), index_t.isin(index_o)
    coverage.update({"events_TAMU": len(ours), "events_gpizzati": len(theirs), "events_in_both_skims": int(in_t.sum())})
    # the other side's row for each event, to compare values and flags
    row_t = pd.Series(np.arange(len(theirs)), index=index_t).reindex(index_o).values
    row_o = pd.Series(np.arange(len(ours)), index=index_o).reindex(index_t).values

    out = args.output_dir / ds
    out.mkdir(parents=True, exist_ok=True)
    plots_dir = args.plots_dir or out / "plots"
    regions = [r for r in REGIONS if not (is_data and r == "Signal_Fit" and not args.unblind)]
    wsum = lambda t, m: float(t.weight.values[m].sum())  # noqa: E731
    yields, cutflows, agreement, lists, factors = [], [], [], [], []
    for category in CATEGORIES:
        for region in regions:
            steps_o, steps_t = passes(ours, category, region), passes(theirs, category, region)
            ff_o, ff_t = first_failed(steps_o), first_failed(steps_t)
            sel_o, sel_t = ff_o == "PASS", ff_t == "PASS"
            both = sel_o & in_t
            both[in_t] &= sel_t[row_t[in_t].astype(int)]
            only_o, only_t = sel_o & ~both, sel_t.copy()
            only_t[in_o] &= ~sel_o[row_o[in_o].astype(int)]
            lost_o = np.where(in_t, ff_t[np.nan_to_num(row_t, nan=0).astype(int)], "not in gpizzati skim")[only_o]
            lost_t = np.where(in_o, ff_o[np.nan_to_num(row_o, nan=0).astype(int)], "not in TAMU skim")[only_t]
            yields.append(dict(category=category, region=region, TAMU=wsum(ours, sel_o), gpizzati=wsum(theirs, sel_t),
                               common=wsum(ours, both), only_TAMU=wsum(ours, only_o), only_gpizzati=wsum(theirs, only_t),
                               # the same numbers restricted to the events present in both skims
                               TAMU_shared=wsum(ours, sel_o & in_t), gpizzati_shared=wsum(theirs, sel_t & in_o),
                               only_TAMU_shared=wsum(ours, only_o & in_t), only_gpizzati_shared=wsum(theirs, only_t & in_o),
                               only_TAMU_lost_by_gpizzati_at=dict(pd.Series(lost_o).value_counts()),
                               only_gpizzati_lost_by_TAMU_at=dict(pd.Series(lost_t).value_counts())))
            # cumulative cutflow on both sides, and events passing on both
            cum_o, cum_t = np.ones(len(ours), bool), np.ones(len(theirs), bool)
            cutflows.append(dict(category=category, region=region, step="skim", TAMU=wsum(ours, cum_o),
                                 gpizzati=wsum(theirs, cum_t), common=wsum(ours, in_t)))
            for step in steps_o:
                cum_o &= steps_o[step]
                cum_t &= steps_t[step]
                common = cum_o & in_t
                common[in_t] &= cum_t[row_t[in_t].astype(int)]
                cutflows.append(dict(category=category, region=region, step=step, TAMU=wsum(ours, cum_o),
                                     gpizzati=wsum(theirs, cum_t), common=wsum(ours, common)))
            # agreement of the variables on the common selected events
            ro, rt = np.where(both)[0], row_t[both].astype(int)
            for name, tolerance in AGREEMENT:
                a, b = ours[name].values[ro], theirs[name].values[rt]
                valid = np.isfinite(a) & np.isfinite(b) & ((b != 0) | (tolerance == 0))
                same = (a == b) if tolerance == 0 else (np.abs(a - b) <= tolerance * np.abs(b))
                relative = (a - b)[valid] / b[valid] if tolerance and valid.any() else np.array([np.nan])
                agreement.append(dict(category=category, region=region, variable=name, tolerance=tolerance,
                                      events=int(valid.sum()), same_fraction=float(same[valid].mean()) if valid.any() else np.nan,
                                      median_rel_diff=float(np.median(relative))))
            if weighted == "analysis" and both.any():
                # weight factors on the common selected events: ratio of the sums and per-event median ratio
                for name in (*FACTORS, "weight"):
                    a, b = ours[name].values[ro], theirs[name].values[rt]
                    with np.errstate(divide="ignore", invalid="ignore"):
                        per_event = a / b
                    factors.append(dict(category=category, region=region, factor=name, events=len(ro),
                                        mean_TAMU=float(np.mean(a)), mean_gpizzati=float(np.mean(b)),
                                        ratio_of_means=float(np.mean(a) / np.mean(b)) if np.mean(b) else np.nan,
                                        median_event_ratio=float(np.nanmedian(per_event)),
                                        fraction_within_1pc=float(np.mean(np.abs(per_event - 1) < 0.01))))
            lists.append(one_side_table(ours, theirs, only_o, row_t, lost_o, "only_TAMU", category, region))
            lists.append(one_side_table(theirs, ours, only_t, row_o, lost_t, "only_gpizzati", category, region))
            if args.plots:
                plot(ours[sel_o], theirs[sel_t], only_o[sel_o], only_t[sel_t], category, region,
                     (args.label_ours, args.label_theirs), plots_dir / f"{category}_{region}")
    table = pd.DataFrame(yields)
    table.drop(columns=["only_TAMU_lost_by_gpizzati_at", "only_gpizzati_lost_by_TAMU_at"]).to_csv(out / "yields.csv", index=False)
    cutflow, agree = pd.DataFrame(cutflows), pd.DataFrame(agreement)
    cutflow.to_csv(out / "cutflow.csv", index=False)
    agree.to_csv(out / "agreement.csv", index=False)
    pd.concat(lists).to_csv(out / "events_one_side_only.csv", index=False)
    factor_table = pd.DataFrame(factors)
    if len(factor_table):
        factor_table.to_csv(out / "weight_factors.csv", index=False)
    (out / "coverage.json").write_text(json.dumps(coverage, indent=2) + "\n")
    write_readme(ds, out, table, cutflow, agree, coverage, args, weighted, plots_dir, factor_table)
    return table, coverage


def one_side_table(this, other, mask, rows, lost, side, category, region):
    """One-side events with the variables and step flags of both skims."""
    columns = ["m_mumu", "mu1_pt", "mu2_pt", "mu1_eta", "mu2_eta", "n_jets", "j1_pt", "j2_pt", "j1_eta", "j2_eta", "mjj"] + \
              [c for c in this.columns if c.startswith("pass_")]
    names = ("TAMU", "gpizzati") if side == "only_TAMU" else ("gpizzati", "TAMU")
    result = this.loc[mask, KEYS].reset_index(drop=True)
    result.insert(0, "side", side)
    result.insert(0, "region", region)
    result.insert(0, "category", category)
    result["lost_by_other_at"] = lost
    mine = this.loc[mask, columns].reset_index(drop=True).add_suffix(f"_{names[0]}")
    idx = rows[mask]
    theirs = other.reindex(np.where(np.isfinite(idx), idx, -1)).loc[:, columns].reset_index(drop=True).add_suffix(f"_{names[1]}")
    return pd.concat([result, mine, theirs], axis=1)


def write_readme(ds, out, table, cutflow, agree, coverage, args, weighted, plots_dir, factor_table=None):
    fmt = {"analysis": lambda x: f"{x:,.4g}", "genWeight": lambda x: f"{x:,.1f}"}.get(weighted, lambda x: f"{x:,.0f}")
    c = coverage
    lines = [f"# Sync {ds} (Run3_2024): TAMU skim_v4 vs gpizzati v7", "",
             f"- TAMU input: `{c['TAMU_input']}` ({c['TAMU_files']} files)",
             f"- gpizzati input: `{c['gpizzati_input']}` ({c['gpizzati_files']} files)",
             f"- {'data' if c['is_data'] else 'MC'}, weights: {weighted or 'none (event counts)'}"
             + (f", files per side limited to {args.max_files}" if args.max_files else ""), "",
             "## Input events (from the skim JSON files)", "",
             "| | JSON found | before skim | after the MC lumi split | after skim |", "|---|---|---|---|---|",
             *[f"| {side} | {t['reports']} | {t['before_skim']:,} | {t['after_lumi_split']:,} | {t['after_skim']:,} |"
               for side, t in c["input_totals"].items()],
             "", "gpizzati's JSON already counts the events after his MC lumi split. TAMU's last skim step is",
             "the number of events written in the skim.", "",
             "## Coverage", "",
             f"- lumisections in both skims: {c['lumis_both']:,}; only TAMU: {c['lumis_only_TAMU']:,}; only gpizzati: {c['lumis_only_gpizzati']:,}",
             f"- events outside the shared lumisections (not compared): TAMU {c['events_TAMU_outside_shared_lumis']:,}, "
             f"gpizzati {c['events_gpizzati_outside_shared_lumis']:,}",
             f"- events with a run:lumi:event used twice (MC extensions), dropped on both sides: "
             f"TAMU {c['events_TAMU_dropped_duplicate_keys']:,}, gpizzati {c['events_gpizzati_dropped_duplicate_keys']:,}",
             f"- skim events compared: TAMU {c['events_TAMU']:,}, gpizzati {c['events_gpizzati']:,}, in both skims {c['events_in_both_skims']:,}", "",
             "## Yields", "",
             "| category | region | TAMU | gpizzati | TAMU / gpizzati | common | only TAMU | only gpizzati |",
             "|---|---|---|---|---|---|---|---|"]
    for r in table.itertuples():
        ratio = f"{r.TAMU / r.gpizzati:.4f}" if r.gpizzati else "-"
        lines.append(f"| {r.category} | {r.region} | {fmt(r.TAMU)} | {fmt(r.gpizzati)} | {ratio} | {fmt(r.common)} | "
                     f"{fmt(r.only_TAMU)} | {fmt(r.only_gpizzati)} |")
    shared = c["events_in_both_skims"] / max(min(c["events_TAMU"], c["events_gpizzati"]), 1)
    lines += ["", "## Yields on the events present in both skims", "",
              "This is the selection comparison proper: the same input events, selected by the two frameworks.", ""]
    if shared < 0.9:
        lines += [f"**Only {shared:.1%} of the skim events are in both skims**: the two productions do not start from the",
                  "same events (for MC, typically a different set of extensions or a different lumi split), so the table",
                  "above compares two different samples and only this one is an event-level sync.", ""]
    lines += ["| category | region | TAMU | gpizzati | TAMU / gpizzati | only TAMU | only gpizzati |", "|---|---|---|---|---|---|---|"]
    for r in table.itertuples():
        ratio = f"{r.TAMU_shared / r.gpizzati_shared:.4f}" if r.gpizzati_shared else "-"
        lines.append(f"| {r.category} | {r.region} | {fmt(r.TAMU_shared)} | {fmt(r.gpizzati_shared)} | {ratio} | "
                     f"{fmt(r.only_TAMU_shared)} | {fmt(r.only_gpizzati_shared)} |")
    lines += ["", "## Where the other side loses the one-side events (event counts)", "",
              "| category | region | only TAMU: gpizzati fails at | only gpizzati: TAMU fails at |", "|---|---|---|---|"]
    for r in table.itertuples():
        show = lambda d: ", ".join(f"{k} {v}" for k, v in d.items()) or "-"  # noqa: E731
        lines.append(f"| {r.category} | {r.region} | {show(r.only_TAMU_lost_by_gpizzati_at)} | {show(r.only_gpizzati_lost_by_TAMU_at)} |")
    for category, region in (("baseline", "Z_sideband"), ("VBF", "Z_sideband")):
        part = cutflow[(cutflow.category == category) & (cutflow.region == region)]
        lines += ["", f"## Cutflow {category}, {region} (cumulative)", "",
                  "| step | TAMU | gpizzati | TAMU / gpizzati | common |", "|---|---|---|---|---|"]
        for r in part.itertuples():
            lines.append(f"| {r.step} | {fmt(r.TAMU)} | {fmt(r.gpizzati)} | {r.TAMU / r.gpizzati:.4f} | {fmt(r.common)} |"
                         if r.gpizzati else f"| {r.step} | {fmt(r.TAMU)} | {fmt(r.gpizzati)} | - | {fmt(r.common)} |")
    if factor_table is not None and len(factor_table):
        lines += ["", "## Analysis weight factors on the common selected events", "",
                  "TAMU / gpizzati. gen = genWeight, pu = pileup, idiso = muon ID x iso SFs, trg = trigger SF,",
                  "norm = lumi x xsec / sum(genWeight) (TAMU: x DY reweights for DY), weight = their product.",
                  f"TAMU events without an analysis weight (not found in the framework output): "
                  f"{c['TAMU_events_without_analysis_weight']:,}.", "",
                  "| category | region | factor | events | mean TAMU | mean gpizzati | ratio of means | median event ratio | events within 1% |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for r in factor_table.itertuples():
            if (r.category, r.region) not in (("baseline", "Z_sideband"), ("baseline", "Signal_Fit"), ("VBF", "Signal_Fit"), ("VBF", "Z_sideband")):
                continue
            lines.append(f"| {r.category} | {r.region} | {r.factor} | {r.events:,} | {r.mean_TAMU:.4g} | {r.mean_gpizzati:.4g} | "
                         f"{r.ratio_of_means:.4f} | {r.median_event_ratio:.4f} | {r.fraction_within_1pc:.3f} |")
    part = agree[(agree.category == "baseline") & (agree.region == "Z_sideband")]
    lines += ["", "## Agreement on the common events (baseline, Z_sideband)", "",
              "| variable | relative tolerance | events | same value | median relative difference |", "|---|---|---|---|---|"]
    for r in part.itertuples():
        lines.append(f"| {r.variable} | {r.tolerance:g} | {r.events:,} | {r.same_fraction:.4f} | "
                     + (f"{r.median_rel_diff:.2e} |" if np.isfinite(r.median_rel_diff) else "exact comparison |"))
    lines += ["", "The trigger step is `muons_trigger` for TAMU and always true for gpizzati (no trigger matching stored).",
              "All categories and regions are in `cutflow.csv` and `agreement.csv`; every one-side event, with the",
              "variables and step flags of both skims, is in `events_one_side_only.csv`."]
    if args.plots:
        lines.append(f"Plots: `{plots_dir}/<category>_<region>/`.")
    (out / "README.md").write_text("\n".join(lines) + "\n")


def main():
    args = arguments()
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})
    for ds in args.dataset:
        table, coverage = compare(ds, args)
        print(f"== {ds}: lumis both {coverage['lumis_both']}, only TAMU {coverage['lumis_only_TAMU']}, "
              f"only gpizzati {coverage['lumis_only_gpizzati']}")
        print(table.drop(columns=["only_TAMU_lost_by_gpizzati_at", "only_gpizzati_lost_by_TAMU_at"]).to_string(index=False))
        print(f"   written {args.output_dir / ds}/README.md", flush=True)


if __name__ == "__main__":
    main()
