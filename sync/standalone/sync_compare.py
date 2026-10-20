#!/usr/bin/env python3
"""Standalone tool to compare the Run 3 H->mumu synchronisation ntuples of several groups.

It needs only numpy, awkward, uproot and matplotlib (see requirements.txt); nothing else from the analysis framework.

Two steps (full instructions in HOWTO.txt):

  1. yields : read the ntuples (events with the same run, lumi, event are counted once) of every group and write tables
       yields_<sample>.csv   events and sum of weights per group, mass region and channel
       overlap_<sample>.csv  events present in both groups of every pair (matched on run, lumi, event)
       mass_<sample>.csv     m(mumu) distribution of every group (selected ggH or VBF events, 70-150 GeV)
  2. plots  : draw bar charts of the yields, the ratio to a reference group, and the overlap matrices from those tables.

  python3 sync_compare.py yields --base /eos/cms/store/group/phys_higgs/cmshmm/run3_sync --outdir out
  python3 sync_compare.py plots  --outdir out --reference TAMU

A group is one entry of GROUPS below: the files of every sample, the tree name, the dimuon-mass branch used for the
regions, one boolean expression per channel (ggH, VBF) and the weight expressions.  A new group is added by adding one entry.
Expressions are Python expressions on branch names; `name[:,0]` is the first entry of a jagged branch.
"""
import argparse
import csv
import glob
import re
import sys
from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

# ----------------------------------------------------------------------------------------------
# Description of the groups (update_after_Hmm_15Oct tuples; paths are relative to --base)
# ----------------------------------------------------------------------------------------------
GROUPS = {'TAMU': {'tree': 'Events',
          'lumi': 'luminosityBlock',
          'mass': 'm_mumu_nano',
          'files': {'Data_Run2024H': ['TAMU/update_after_Hmm_15Oct/TAMU-Run2024H.root'],
                    'ggH': ['TAMU/update_after_Hmm_15Oct/TAMU-GluGluHto2Mu.root'],
                    'VBFH': ['TAMU/update_after_Hmm_15Oct/TAMU-VBFHto2Mu_M125_powheg.root',
                             'TAMU/update_after_Hmm_15Oct/TAMU-VBFHto2Mu_M125_powheg_ext1.root']},
          'channels': {'ggH': 'ggF', 'VBF': 'VBF'},
          'weights': {'gen': 'weight_generator', 'all': 'weight_generator*weight_pileup*weight_muon_idiso*weight_trigger'}},
 'Pisa': {'tree': 'Events',
          'lumi': 'luminosityBlock',
          'mass': 'Dimuon_rawMass',
          'files': {'Data_Run2024H': ['PISA/update_after_Hmm_15Oct/PISA-Muon0_Run2024H-MINIv6NANOv15-v1/Events.root',
                                      'PISA/update_after_Hmm_15Oct/PISA-Muon1_Run2024H-MINIv6NANOv15-v2/Events.root'],
                    'ggH': ['PISA/update_after_Hmm_15Oct/PISA-GluGluH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15/Events.root'],
                    'VBFH': ['PISA/update_after_Hmm_15Oct/PISA-VBFH-Hto2Mu_Par-M-125_2024Summer24NanoAODv15/Events.root']},
          'channels': {'VBF': 'pass_BJetVeto', 'ggH': 'pass_EleVeto & ~pass_DijetDeltaEta & (nMediumBTagJet == 0) & (nLooseBTagJet < 2)'},
          'weights': {'gen': 'Generator_weight', 'all': 'Generator_weight*PuSF*MuIdSF*MuIsoSF*MuTrgSF'}},
 'Marmara': {'tree': 'HmmTree',
             'lumi': 'lumi',
             'mass': 'mm_mass_raw',
             'files': {'Data_Run2024H': ['Marmara/update_after_Hmm_15Oct/tree_DATA_2024_loose_round2_update_last.root'],
                       'ggH': ['Marmara/update_after_Hmm_15Oct/tree_GGH_2024_loose_round2_upd2_nobtag_weight.root'],
                       'VBFH': ['Marmara/update_after_Hmm_15Oct/tree_VBF_2024_loose_round2_upd2_nobtag_weight.root']},
             'channels': {'ggH': 'sel_full & cat_isGGH[:,0]', 'VBF': 'sel_full & cat_isVBF[:,0]'},
             'weights': {'gen': 'w_gen', 'all': 'w_gen*w_norm*w_pu*w_muID_tot'}},
 'Caltech': {'tree': 'tree_output',
             'lumi': 'luminosityBlock',
             'mass': 'diMuon_bscMass_noCorr',
             'files': {'Data_Run2024H': ['Caltech/update_after_Hmm_15Oct/Muon0_2024H_skim.root',
                                         'Caltech/update_after_Hmm_15Oct/Muon1_2024H_skim.root'],
                       'ggH': ['Caltech/update_after_Hmm_15Oct/ggH_skim.root'],
                       'VBFH': ['Caltech/update_after_Hmm_15Oct/VBF_skim.root']},
             'channels': {'ggH': '(is_ggH_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto',
                          'VBF': '(is_VBF_category == 1) & event_pass_muon_trigger_match & event_pass_jet_veto'},
             'weights': {'gen': 'gen_weight', 'all': 'gen_weight*pileup_weight*mu_trig_eff_SF*mu_ID_eff_SF*mu_ISO_eff_SF'}},
 'NNU': {'tree': 'Events',
         'lumi': 'luminosityBlock',
         'mass': 'dimuon_raw_mass',
         'files': {'Data_Run2024H': ['NNU/update_after_Hmm_15Oct/NNU-dataH2024.root'],
                   'VBFH': ['NNU/update_after_Hmm_15Oct/NNU-VBFH-Hto2Mu_M125_RunIII2024Summer24_20261020.root']},
         'channels': {'VBF': 'passAnalysisSelection',
                      'ggH': '(passLumiJSON & passTrigger & passPV & passMETFilters & passTriggerMatch & passRawPreselection & '
                             'passMuonSelection & passThirdLeptonVeto & passBTagVeto) & (passAnalysisSelection==0)'},
         'weights': {'gen': 'genWeight', 'all': 'weightBaseNoNorm'}},
 'Iowa': {'tree': 'SyncEvents',
          'lumi': 'luminosityBlock',
          'mass': 'massRaw',
          'files': {'Data_Run2024H': ['Iowa/update_after_Hmm_15Oct/Muon0/*.root', 'Iowa/update_after_Hmm_15Oct/Muon1/*.root'],
                    'VBFH': ['Iowa/update_after_Hmm_15Oct/VBF/*.root']},
          'channels': {'ggH': 'passGolden & passFilters & passHLTAndObjectMatch & passAnalysisMuonSelection & passExactlyTwoAnalysisMuons '
                              '& passElectronVeto & passAdditionalLeptonVeto & passJetVetoMap & passBJetVeto & passVertex & '
                              'analysisPairPassesRawWindow & oppositeSign & ~passVBFCategoryApprox',
                       'VBF': 'passGolden & passFilters & passHLTAndObjectMatch & passAnalysisMuonSelection & passExactlyTwoAnalysisMuons '
                              '& passElectronVeto & passAdditionalLeptonVeto & passJetVetoMap & passBJetVeto & passVertex & '
                              'analysisPairPassesRawWindow & oppositeSign & passVBFCategoryApprox'},
          'weights': {'gen': 'weightGenerator', 'all': 'weightGenerator*weightPileup*weightMuonIDISO*weightTrigger'}}}

SAMPLES = ["Data_Run2024H", "ggH", "VBFH"]
CHANNELS = ["ggH", "VBF"]
# mass regions on the group's raw dimuon mass; the data are blinded in the signal window
REGIONS = {"Z_sideband": (70, 110), "H_sideband": None, "Signal_Fit": (115, 135)}
H_SIDEBAND = ((110, 115), (135, 150))
MASS_RANGE = (70, 150)


def in_region(mass, region):
    if region == "H_sideband":
        return ((mass > 110) & (mass < 115)) | ((mass > 135) & (mass < 150))
    low, high = REGIONS[region]
    return (mass > low) & (mass < high)


# ----------------------------------------------------------------------------------------------
# Reading the ntuples
# ----------------------------------------------------------------------------------------------

def first(array):
    """First entry of a jagged branch, False/NaN when missing."""
    padded = ak.pad_none(array, 1, axis=1)[:, 0]
    return ak.fill_none(padded, False if "bool" in str(ak.type(array)) else np.nan)


def evaluate(arrays, expression):
    """Evaluate an expression on the branches of a chunk and return a flat numpy array."""
    expression = re.sub(r"(\w+)\[:,0\]", r"first(\1)", expression)
    scope = {name: arrays[name] for name in branch_names(expression) if name in arrays.fields}
    scope["first"] = first
    return np.asarray(ak.to_numpy(eval(expression, {"__builtins__": {}}, scope)))


def branch_names(expression):
    expression = re.sub(r"(\w+)\[:,0\]", r"first(\1)", expression)
    return [n for n in dict.fromkeys(re.findall(r"[A-Za-z_]\w*", expression)) if n != "first"]


def sample_files(base, spec, sample):
    paths = []
    for pattern in spec["files"][sample]:
        paths += sorted(glob.glob(str(Path(base) / pattern)))
    if not paths:
        raise FileNotFoundError(f"no file for {sample}: {spec['files'][sample]} under {base}")
    return paths


def pack_keys(run, lumi, event, is_data):
    """One integer per event (run, lumi, event), used to match events between groups."""
    run, lumi, event = (np.asarray(x).astype(np.uint64) for x in (run, lumi, event))
    if is_data:
        return (run - np.uint64(370000)) << np.uint64(48) | lumi << np.uint64(34) | event
    return np.uint64(1) << np.uint64(63) | run << np.uint64(59) | lumi << np.uint64(39) | event


def read_sample(base, spec, sample, step_size=500_000):
    """Selected events of one group and sample: key, mass, channel flags and weights (70 < m < 150 GeV only)."""
    is_data = sample.startswith("Data")
    expressions = [spec["mass"], *spec["channels"].values()] + ([] if is_data else list(spec["weights"].values()))
    wanted = ["run", spec["lumi"], "event"]
    for expression in expressions:
        wanted += branch_names(expression)
    out = {k: [] for k in ("key", "mass", "ggH", "VBF", "w_gen", "w_all")}
    entries = 0
    for path in sample_files(base, spec, sample):
        tree = uproot.open(path)[spec["tree"]]
        names = [n for n in dict.fromkeys(wanted) if n in tree.keys()]
        for chunk in tree.iterate(names, step_size=step_size, library="ak"):
            entries += len(chunk)
            mass = evaluate(chunk, spec["mass"]).astype(np.float64)
            keep = (mass > MASS_RANGE[0]) & (mass < MASS_RANGE[1])
            channel = {c: evaluate(chunk, spec["channels"][c]).astype(bool) for c in CHANNELS}
            keep &= channel["ggH"] | channel["VBF"]
            out["key"].append(pack_keys(chunk["run"], chunk[spec["lumi"]], chunk["event"], is_data)[keep])
            out["mass"].append(mass[keep])
            for c in CHANNELS:
                out[c].append(channel[c][keep])
            for mode, name in (("w_gen", "gen"), ("w_all", "all")):
                w = np.ones(len(chunk)) if is_data else evaluate(chunk, spec["weights"][name]).astype(np.float64)
                out[mode].append(w[keep])
    result = {k: np.concatenate(v) for k, v in out.items()}
    result["entries"] = np.array(entries)
    return result


def load_or_read(base, group, sample, outdir):
    """Read the ntuple once; the selected events are cached in <outdir>/cache/<group>_<sample>.npz."""
    cache = Path(outdir) / "cache" / f"{group}_{sample}.npz"
    if cache.exists():
        return dict(np.load(cache))
    print(f"[read] {group} {sample}", flush=True)
    result = read_sample(base, GROUPS[group], sample)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **result)
    return result


def deduplicate(d):
    """Count every (run, lumi, event) once per group and sample (some files contain events twice).
    `duplicates` is the number of rows dropped; for VBFH (main + ext1) a few identical identifiers of distinct events (<0.01%) are dropped too."""
    _, index = np.unique(d["key"], return_index=True)
    index = np.sort(index)
    out = {k: v[index] for k, v in d.items() if k != "entries"}
    out["entries"] = d["entries"]
    out["duplicates"] = np.array(len(d["key"]) - len(index))
    return out


# ----------------------------------------------------------------------------------------------
# Step 1: tables
# ----------------------------------------------------------------------------------------------

def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def cells(sample):
    """(region, channel) pairs of a sample; the data signal window is blinded."""
    return [(r, c) for r in REGIONS for c in CHANNELS if not (sample.startswith("Data") and r == "Signal_Fit")]


def step_yields(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    groups = args.groups or list(GROUPS)
    for sample in args.samples or SAMPLES:
        data = {g: deduplicate(load_or_read(args.base, g, sample, outdir)) for g in groups if sample in GROUPS[g]["files"]}
        yields, overlap, mass_rows = [], [], []
        for g, d in data.items():
            for region, channel in cells(sample):
                sel = d[channel] & in_region(d["mass"], region)
                yields.append([g, sample, region, channel, int(sel.sum()), f"{d['w_gen'][sel].sum():.6g}", f"{d['w_all'][sel].sum():.6g}", int(d["duplicates"])])
            for channel in CHANNELS:
                sel = d[channel]
                counts, edges = np.histogram(d["mass"][sel], bins=np.arange(MASS_RANGE[0], MASS_RANGE[1] + 1, 1.0))
                mass_rows += [[g, channel, edges[i], edges[i + 1], int(counts[i])] for i in range(len(counts))]
        names = list(data)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                for region, channel in cells(sample):
                    ka = np.unique(data[a]["key"][data[a][channel] & in_region(data[a]["mass"], region)])
                    kb = np.unique(data[b]["key"][data[b][channel] & in_region(data[b]["mass"], region)])
                    common = len(np.intersect1d(ka, kb, assume_unique=True))
                    overlap.append([a, b, sample, region, channel, len(ka), len(kb), common,
                                    f"{common / max(len(ka), 1):.5f}", f"{common / max(len(kb), 1):.5f}"])
        write_csv(outdir / f"yields_{sample}.csv", ["group", "sample", "region", "channel", "events", "sum_gen_weight", "sum_all_weights", "duplicate_rows_dropped"], yields)
        write_csv(outdir / f"overlap_{sample}.csv", ["groupA", "groupB", "sample", "region", "channel", "nA", "nB", "common", "common_over_A", "common_over_B"], overlap)
        write_csv(outdir / f"mass_{sample}.csv", ["group", "channel", "low", "high", "events"], mass_rows)
        print(f"[tables] {sample}: {len(yields)} yield rows, {len(overlap)} overlap rows", flush=True)


# ----------------------------------------------------------------------------------------------
# Step 2: plots
# ----------------------------------------------------------------------------------------------

def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def colours(groups):
    return {g: plt.cm.tab10(i % 10) for i, g in enumerate(groups)}


def plot_bars(outdir, sample, rows):
    groups = list(dict.fromkeys(r["group"] for r in rows))
    labels = [f"{r}\n{c}" for r, c in cells(sample)]
    colour = colours(groups)
    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.8 / len(groups)
    for i, g in enumerate(groups):
        y = [next(int(r["events"]) for r in rows if (r["group"], r["region"], r["channel"]) == (g, reg, ch)) for reg, ch in cells(sample)]
        ax.bar(np.arange(len(labels)) + (i - (len(groups) - 1) / 2) * width, y, width, label=g, color=colour[g])
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("selected events")
    ax.set_title(f"{sample}: events per mass region and channel")
    ax.legend(ncol=len(groups), fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"yield_bars_{sample}.png", dpi=130)
    plt.close(fig)


def plot_ratio(outdir, sample, rows, reference):
    groups = [g for g in dict.fromkeys(r["group"] for r in rows) if g != reference]
    if reference not in {r["group"] for r in rows}:
        return
    labels = [f"{r}\n{c}" for r, c in cells(sample)]
    colour = colours([reference] + groups)
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.axhspan(0.99, 1.01, color="0.85", label="±1%")
    ax.axhspan(0.95, 1.05, color="0.93", zorder=0, label="±5%")
    for i, g in enumerate(groups):
        ratios = []
        for reg, ch in cells(sample):
            number = {r["group"]: int(r["events"]) for r in rows if (r["region"], r["channel"]) == (reg, ch)}
            ratios.append(number[g] / number[reference] if number.get(reference) else np.nan)
        x = np.arange(len(labels)) + (i - (len(groups) - 1) / 2) * 0.12
        ax.plot(x, ratios, "o", color=colour[g], label=g)
    ax.axhline(1, color="k", lw=0.8)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 3)
    ax.set_yticks([0.1, 0.2, 0.5, 0.8, 0.9, 1, 1.1, 1.25, 2, 3])
    ax.set_yticklabels(["0.1", "0.2", "0.5", "0.8", "0.9", "1", "1.1", "1.25", "2", "3"])
    ax.set_ylabel(f"events / {reference}")
    ax.set_title(f"{sample}: yield agreement with {reference}")
    ax.legend(ncol=len(groups) + 2, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"yield_ratio_{sample}.png", dpi=130)
    plt.close(fig)


def plot_overlap(outdir, sample, rows):
    groups = sorted({r["groupA"] for r in rows} | {r["groupB"] for r in rows})
    if len(groups) < 2:
        return
    fig, axes = plt.subplots(1, len(CHANNELS), figsize=(6 * len(CHANNELS), 5.2))
    for ax, channel in zip(axes, CHANNELS):
        matrix = np.full((len(groups), len(groups)), np.nan)
        for a_i, a in enumerate(groups):
            for b_i, b in enumerate(groups):
                if a == b:
                    continue
                nA = nB = common = 0
                for r in rows:
                    if r["channel"] == channel and {r["groupA"], r["groupB"]} == {a, b}:
                        nA += int(r["nA"]) if r["groupA"] == a else int(r["nB"])
                        nB += int(r["nB"]) if r["groupA"] == a else int(r["nA"])
                        common += int(r["common"])
                if nA and nB:
                    matrix[a_i, b_i] = common / min(nA, nB)
        image = ax.imshow(matrix, vmin=0.7, vmax=1.0, cmap="viridis")
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=45, ha="right")
        ax.set_yticks(range(len(groups)))
        ax.set_yticklabels(groups)
        for i in range(len(groups)):
            for j in range(len(groups)):
                if not np.isnan(matrix[i, j]):
                    ax.text(j, i, f"{100 * matrix[i, j]:.1f}", ha="center", va="center", color="w", fontsize=8)
        ax.set_title(f"{sample} {channel}: common / smaller selection [%]")
        fig.colorbar(image, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(outdir / f"overlap_matrix_{sample}.png", dpi=130)
    plt.close(fig)


def plot_mass(outdir, sample, rows):
    fig, axes = plt.subplots(1, len(CHANNELS), figsize=(6 * len(CHANNELS), 4.5))
    groups = list(dict.fromkeys(r["group"] for r in rows))
    colour = colours(groups)
    for ax, channel in zip(axes, CHANNELS):
        for g in groups:
            sel = [r for r in rows if r["group"] == g and r["channel"] == channel]
            if not sel:
                continue
            edges = np.array([float(r["low"]) for r in sel] + [float(sel[-1]["high"])])
            counts = np.array([float(r["events"]) for r in sel])
            ax.stairs(counts, edges, label=g, color=colour[g])
        ax.set_xlabel("raw m(mumu) [GeV]")
        ax.set_ylabel("events / GeV")
        ax.set_title(f"{sample} {channel}")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"mass_{sample}.png", dpi=130)
    plt.close(fig)


def step_plots(args):
    outdir = Path(args.outdir)
    for sample in args.samples or SAMPLES:
        if (outdir / f"yields_{sample}.csv").exists():
            rows = read_csv(outdir / f"yields_{sample}.csv")
            plot_bars(outdir, sample, rows)
            plot_ratio(outdir, sample, rows, args.reference)
            plot_overlap(outdir, sample, read_csv(outdir / f"overlap_{sample}.csv"))
            plot_mass(outdir, sample, read_csv(outdir / f"mass_{sample}.csv"))
            print(f"[plots] {sample}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="step", required=True)
    for name in ("yields", "plots"):
        p = sub.add_parser(name)
        p.add_argument("--outdir", default="sync_out", help="where tables, plots and the cache are written")
        p.add_argument("--samples", nargs="+", choices=SAMPLES)
        if name == "yields":
            p.add_argument("--base", default="/eos/cms/store/group/phys_higgs/cmshmm/run3_sync", help="directory with one folder per group")
            p.add_argument("--groups", nargs="+", choices=list(GROUPS))
        else:
            p.add_argument("--reference", default="TAMU", help="group the yield ratios are computed against")
    args = parser.parse_args()
    {"yields": step_yields, "plots": step_plots}[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
