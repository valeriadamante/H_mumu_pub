#!/usr/bin/env python3
"""Collect the correction configuration actually wired into the H->mumu skim.

The inventory is built statically from the repository sources: importing the
correction modules would require ROOT, correctionlib and valid CVMFS payloads.
By default a human-readable Markdown report is printed; ``--format json`` is
available for machine-readable output.
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path
from typing import Any

import yaml


def literal_assignments(path: Path) -> dict[str, Any]:
    """Return every assignment whose value can be evaluated as a literal."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    values: dict[str, Any] = {}
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        try:
            value = ast.literal_eval(node.value)
        except (ValueError, TypeError):
            continue
        for target in targets:
            if isinstance(target, ast.Name):
                values[target.id] = value
    return values


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        value = yaml.safe_load(stream) or {}
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a top-level mapping")
    return value


def cvmfs_path(template: str, folder: str, **fields: str) -> str:
    return template.format(folder, **fields)


def build_inventory(repo: Path, requested_eras: list[str]) -> dict[str, Any]:
    corr = repo / "corrections"
    general = literal_assignments(corr / "general.py")
    jets = literal_assignments(corr / "jets.py")
    pu = literal_assignments(corr / "pu.py")
    mu = literal_assignments(corr / "mu.py")
    veto = literal_assignments(corr / "jetVetoMap.py")
    btag = literal_assignments(corr / "btag_wpValues.py")

    period_names = general["period_names"]
    pog = general["pog_folder_names"]
    known_eras = sorted(
        p.name
        for p in (repo / "config").glob("Run3_*")
        if p.is_dir() and p.name in period_names
    )
    eras = requested_eras or known_eras
    unknown = sorted(set(eras) - set(period_names))
    if unknown:
        raise ValueError(f"unknown era(s): {', '.join(unknown)}")

    pu_template = "/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/{folder}/latest/puWeights{suffix}.json.gz"
    mu_template = "/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/{}/latest/muon_Z.json.gz"
    veto_template = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{}/latest/jetvetomaps.json.gz"
    scare_template = "corrections/data/MUO/MuonScaRe/{}/{}"

    result: dict[str, Any] = {
        "repository": str(repo),
        "source_files": [
            "analysis/skim.py", "corrections/general.py", "corrections/jets.py",
            "corrections/pu.py", "corrections/mu.py", "corrections/muon_scare.py",
            "corrections/muon_fsr.py", "corrections/jetVetoMap.py",
            "corrections/btag_wpValues.py",
        ],
        "correction_chain": {
            "data": ["Golden JSON", "Muon ScaRe", "Muon FSR", "JEC", "jet veto map", "b-tag working points"],
            "mc": ["pileup", "cross section/luminosity", "Muon ScaRe", "Muon FSR", "JEC", "JER", "muon ID/ISO/trigger SF", "jet veto map", "b-tag working points"],
        },
        "eras": {},
    }

    for era in eras:
        cfg_path = repo / "config" / era / "maincfg.yaml"
        trg_path = repo / "config" / era / "triggers.yaml"
        syst_path = repo / "config" / era / "systematics.yaml"
        cfg = load_yaml(cfg_path)
        trg = load_yaml(trg_path)
        syst = load_yaml(syst_path)
        period = period_names[era]
        jerc_folder = pog["JERC"][period]
        btv_folder = pog["BTV"][period]
        muo_folder = pog["MUO"][period]
        lum_folder = pog["LUM"][period]
        suffix = "_BCDEFGHI" if period == "2024_Summer24" else ""
        if period in {"2025_Summer24", "2026_Summer24"}:
            suffix = "_2025pp_Golden_Summer24_25ns_69200ub"

        trigger_keys = []
        for trigger in trg.values():
            for leg in trigger.get("legs", []) if isinstance(trigger, dict) else []:
                key = leg.get("jsonTRGcorrection_key") if isinstance(leg, dict) else None
                if key:
                    trigger_keys.append(key)

        mu_sources = (
            mu.get("MediumMuIDIso_SF_Sources", {}).get(period, [])
            + mu.get("MediumMuReco_SF_sources", {}).get(period, [])
            + mu.get("MediumMuTrg_SF_Sources", {}).get(period, [])
        )
        requested = cfg.get("requested_SFs", [])
        if requested:
            names = mu["MediumMu_SF_Sources_dict"]
            mu_sources = [key for key in mu_sources if names.get(key) in requested]

        mu_sf_path = mu_template.format(muo_folder)
        if era == "Run3_2025":
            mu_sf_path = "corrections/data/MUO/SF/Run3-25Prompt-Summer24-NanoAODv15/muon_Z.json.gz"

        result["eras"][era] = {
            "period": period,
            "nanoAOD": cfg.get("nano_version"),
            "variations": bool(cfg.get("want_variations", False)),
            "golden_json": cfg.get("lumiFile"),
            "trigger_paths": [path for item in trg.values() if isinstance(item, dict) for path in item.get("path", [])],
            "trigger_sf_keys": trigger_keys,
            "jets": {
                "algorithm": jets["jet_algorithm"],
                "apply_JEC": bool(cfg.get("apply_JES", True)),
                "apply_JER_mc": bool(cfg.get("apply_JER", True)),
                "json": jets["jet_jsonPath"].format(jerc_folder),
                "smearing_json": jets["jetsmear_jsonFile"],
                "JEC_MC_tags": jets["jec_tag_map_mc"].get(period, []),
                "JEC_DATA_tag_templates": jets["jec_tag_map_data"].get(period, []),
                "JER_tag": jets["jer_tag_map"].get(period),
                "uncertainty_sources": (["JER"] + (jets["unc_sources_regrouped"] if cfg.get("use_regrouped", False) else jets["uncSources_minimal"])),
            },
            "pileup": {
                "json": pu_template.format(folder=lum_folder, suffix=suffix),
                "key": pu["golden_json_dict"].get(period),
            },
            "muons": {
                "SF_json": mu_sf_path,
                "SF_keys": mu_sources,
                "ScaRe_json": scare_template.format(muo_folder, "muon_scalesmearing.json"),
                "ScaRe_VXBS_json": scare_template.format(muo_folder, "muon_scalesmearing_VXBS.json"),
                "FSR": "algorithmic correction in corrections/muon_fsr.py (no payload)",
            },
            "jet_veto_map": {
                "json": veto_template.format(jerc_folder),
                "key": veto.get("jetvetomap_names", {}).get(period),
            },
            "btag_working_points": {
                "tagger": cfg.get("bTagAlgo", "PNet"),
                "json": btag["JSON_PATH"].format(btv_folder),
                "fallback_json": btag["ALTERNATIVE_JSON_PATH"].format(btv_folder),
                "key": f"{cfg.get('bTagAlgo', 'PNet')}_wp_values",
                "working_points": ["L", "M", "T"],
            },
            "lhe_theory_weights": {
                "PDF": syst.get("pdf", {}),
                "QCD_scale": syst.get("qcd_scale", {}),
            },
        }
    return result


def markdown(data: dict[str, Any]) -> str:
    lines = ["# Correction inventory", "", "Generated from the analysis configuration and correction source maps.", ""]
    for era, item in data["eras"].items():
        jets = item["jets"]
        lines += [f"## {era}", "", f"- Period: `{item['period']}`; NanoAOD: `{item['nanoAOD']}`; variations: `{item['variations']}`", f"- Golden JSON: `{item['golden_json']}`", f"- Trigger: `{', '.join(item['trigger_paths'])}`; SF key: `{', '.join(item['trigger_sf_keys'])}`", f"- JEC/JER JSON: `{jets['json']}`", f"- JEC MC tags: `{', '.join(jets['JEC_MC_tags'])}`", f"- JEC DATA tag templates: `{', '.join(jets['JEC_DATA_tag_templates'])}`", f"- JER tag: `{jets['JER_tag']}`; algorithm: `{jets['algorithm']}`", f"- Jet uncertainties: `{', '.join(jets['uncertainty_sources'])}`", f"- JER smearing JSON: `{jets['smearing_json']}`", f"- Pileup: `{item['pileup']['json']}`; key: `{item['pileup']['key']}`", f"- Muon SF: `{item['muons']['SF_json']}`", f"- Muon SF keys: `{', '.join(item['muons']['SF_keys'])}`", f"- Muon ScaRe: `{item['muons']['ScaRe_json']}`; VXBS: `{item['muons']['ScaRe_VXBS_json']}`", f"- Jet veto map: `{item['jet_veto_map']['json']}`; key: `{item['jet_veto_map']['key']}`", f"- b-tag WP ({item['btag_working_points']['tagger']}): `{item['btag_working_points']['json']}`; key: `{item['btag_working_points']['key']}`", ""]
    return "\n".join(lines)


def tex_escape(value: Any) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def latex(data: dict[str, Any]) -> str:
    """Create an includable longtable with one row per correction and era."""
    rows: list[tuple[str, str, str, str]] = []
    for era, item in data["eras"].items():
        jets = item["jets"]
        pdf = item["lhe_theory_weights"]["PDF"]
        qcd = item["lhe_theory_weights"]["QCD_scale"]
        replica_last = int(pdf.get("replica_start", 1)) + int(pdf.get("replica_count", 100)) - 1
        qcd_points = "; ".join(
            f"{point['name']}=index {point['index']}"
            for point in qcd.get("points", [])
        )
        qcd_pairs = "; ".join(
            f"{variation['name']}: {variation['down']} / {variation['up']}"
            for variation in qcd.get("variations", [])
        )
        rows.extend(
            [
                (era, "JEC/JER", jets["json"],
                 "JEC MC: " + ", ".join(jets["JEC_MC_tags"])
                 + "; JEC data: " + ", ".join(jets["JEC_DATA_tag_templates"])
                 + "; JER: " + str(jets["JER_tag"])),
                (era, "JER smearing", jets["smearing_json"],
                 "sources: " + ", ".join(jets["uncertainty_sources"])),
                (era, "Pileup", item["pileup"]["json"], item["pileup"]["key"]),
                (era, "Muon SF", item["muons"]["SF_json"],
                 ", ".join(item["muons"]["SF_keys"])),
                (era, "Muon ScaRe", item["muons"]["ScaRe_json"],
                 "standard; VXBS payload: " + item["muons"]["ScaRe_VXBS_json"]),
                (era, "Muon FSR", "---", "algorithmic correction; no external payload"),
                (era, "Jet veto map", item["jet_veto_map"]["json"], item["jet_veto_map"]["key"]),
                (era, "b-tag WP", item["btag_working_points"]["json"],
                 item["btag_working_points"]["key"] + " (L, M, T)"),
                (era, "Golden JSON", item["golden_json"], "data luminosity mask"),
                (era, "Trigger SF", ", ".join(item["trigger_paths"]),
                 ", ".join(item["trigger_sf_keys"])),
                (era, "LHE PDF uncertainty", pdf.get("branch", "LHEPdfWeight"),
                 f"{pdf.get('prescription', 'replicas_rms')}: nominal index "
                 f"{pdf.get('nominal_index', 0)}; replicas "
                 f"{pdf.get('replica_start', 1)}--{replica_last}; "
                 "up=1+RMS, down=max(0,1-RMS); alphaS excluded; "
                 f"missing branch={pdf.get('missing_branch', 'unity')}"),
                (era, "LHE QCD scales", qcd.get("branch", "LHEScaleWeight"),
                 qcd_points + "; dataset-normalized six-point prescription; "
                 + qcd_pairs + "; denominator sums are stored in skim reports; "
                 f"missing sums={qcd.get('missing_sums', 'error')}"),
            ]
        )

    lines = [
        r"% Generated by sync/python/collect_corrections.py; do not edit manually.",
        r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.10\textwidth} >{\raggedright\arraybackslash}p{0.10\textwidth} >{\raggedright\arraybackslash}p{0.38\textwidth} >{\raggedright\arraybackslash}p{0.37\textwidth}}",
        r"\caption{Correction payloads, paths, and tags used in the analysis.}\label{tab:corrections} \\",
        r"\toprule",
        r"\textbf{Era} & \textbf{Correction} & \textbf{Payload / path} & \textbf{Key / tag} \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{4}{c}{\tablename\ \thetable\ -- continued} \\",
        r"\toprule",
        r"\textbf{Era} & \textbf{Correction} & \textbf{Payload / path} & \textbf{Key / tag} \\",
        r"\midrule",
        r"\endhead",
        r"\midrule \multicolumn{4}{r}{Continued on next page} \\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
    ]
    previous_era = None
    for era, correction, payload, key in rows:
        era_cell = tex_escape(era) if era != previous_era else ""
        lines.append(
            f"{era_cell} & {tex_escape(correction)} & "
            rf"\path{{{payload}}} & \texttt{{{tex_escape(key)}}} \\"
        )
        previous_era = era
    lines += [
        r"\end{longtable}",
        "",
        r"\paragraph{LHE PDF uncertainty.}",
        r"The nominal member is \texttt{LHEPdfWeight[0]} and replicas 1--100 are used. "
        r"For each event, the relative replica RMS is "
        r"$\sigma_{\mathrm{PDF}}=\sqrt{\sum_i(w_i/w_0-1)^2/(N_{\mathrm{valid}}-1)}$. "
        r"The multiplicative variations are $w_{\mathrm{PDF}}^{\mathrm{Up}}=1+\sigma_{\mathrm{PDF}}$ "
        r"and $w_{\mathrm{PDF}}^{\mathrm{Down}}=\max(0,1-\sigma_{\mathrm{PDF}})$. "
        r"The configured alphaS members are not included.",
        "",
        r"\paragraph{LHE QCD-scale uncertainty.}",
        r"$\mu_R$ is the QCD renormalization scale, which controls the scale at which "
        r"$\alpha_s$ is evaluated. $\mu_F$ is the factorization scale, which separates "
        r"long-distance PDF evolution from the short-distance hard process. Values 0.5, 1, "
        r"and 2 are multiplicative factors relative to the nominal generator scales. "
        r"The configured final nuisances are: "
        r"\texttt{QCD\_ren\_scale}: Down=$(\mu_R,\mu_F)=(0.5,1)$ "
        r"(\texttt{LHEScaleWeight[1]}), Up=$(2,1)$ (index 7); and "
        r"\texttt{QCD\_fac\_scale}: Down=$(1,0.5)$ (index 3), Up=$(1,2)$ (index 5). "
        r"The corner points $(0.5,0.5)$ (index 0) and $(2,2)$ (index 8) are listed and "
        r"normalized in the skim metadata but are not used by the two final Up/Down pairs. "
        r"Each scale template is normalized with its own dataset-wide denominator from the "
        r"skim reports before histogram filling.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", action="append", default=[], help="era to include (repeatable; default: all Run3 eras)")
    parser.add_argument("--format", choices=("markdown", "json", "latex"), default="markdown")
    parser.add_argument("--output", type=Path, help="write the report to this file instead of stdout")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    try:
        inventory = build_inventory(repo, args.era)
    except (KeyError, ValueError, OSError, SyntaxError) as exc:
        parser.error(str(exc))
    renderers = {"markdown": markdown, "latex": latex}
    text = json.dumps(inventory, indent=2) + "\n" if args.format == "json" else renderers[args.format](inventory) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
