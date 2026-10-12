#!/usr/bin/env python3
"""Generate the LaTeX inventory of selections applied by the NanoAOD skim.

Configuration-dependent values are read from every supported Run3 YAML file.
The execution order and the presence/absence of event filters are checked
against the Python sources, so an incompatible workflow change stops generation
instead of silently producing stale documentation.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import yaml


ERAS = (
    "Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix",
    "Run3_2024", "Run3_2025", "Run3_2026",
)


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        data = yaml.safe_load(stream) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a top-level mapping")
    return data


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def require(source: str, pattern: str, location: str) -> None:
    if re.search(pattern, source, re.MULTILINE | re.DOTALL) is None:
        raise RuntimeError(f"workflow check failed in {location}: {pattern!r}")


def validate_workflow(repo: Path) -> None:
    """Guard the semantic assumptions used to describe the skim."""
    skim = read(repo / "analysis/skim.py")
    ordered_markers = [
        "ApplyOrthogonalLumiFilter", "ApplyGenVBFFilter", "define_base_weights",
        "apply_golden_json", "apply_corrections", "applyMETFlags",
        "ApplyMuonTriggerMatching", "ProcessMuonVariables", "ApplyElectronVeto",
        "apply_muIDIso_weights", "ProcessAllJetVariables", "ApplyJetVetoMap",
        "SelectJetVars", "SelectVBFJets", '.Snapshot("Events"',
    ]
    # rfind skips the import occurrence and selects the actual invocation.
    positions = [skim.rfind(marker) for marker in ordered_markers]
    if any(position < 0 for position in positions) or positions != sorted(positions):
        raise RuntimeError("analysis/skim.py workflow changed: update the table generator")

    muons = read(repo / "analysis/muons.py")
    require(muons, r'good_muons\{suff\}",f"Muon_pt_raw_noCorr > \{pt_min\} && abs\(Muon_eta\) < 2\.4 && Muon_mediumId && Muon_pfIsoId >= 2', "analysis/muons.py")
    # Trigger matching and the muon SFs use the uncorrected pT (rows 9 and 14).
    require(muons, r'pt = "pt_raw_noCorr"\s+muon_p4 = "Muon_p4_raw_noCorr"', "analysis/muons.py")
    require(read(repo / "corrections/mu.py"), r'p4_pt = f"mu\{leg_idx\}_pt_raw_noCorr"', "corrections/mu.py")
    require(muons, r'Filter\(" && "\.join\(event_filters\), "Exactly 2 muons"\)', "analysis/muons.py")
    require(muons, r'Filter\(" && "\.join\(mass_filters\), "dimuon mass cut"\)', "analysis/muons.py")
    require(muons, r'Electron_pt > 20 && abs\(Electron_eta\) < 2\.5 && Electron_mvaIso_WP90', "analysis/muons.py")

    jets = read(repo / "analysis/jets.py")
    require(jets, r'Jet_preSel\{suff\}.*?> \{pt_min\}.*?< \{eta_max\}.*?Jet_passJetIdTight', "analysis/jets.py")
    require(jets, r'DeltaR\(\{p4_branch\}\[i\], mu1_p4\) > 0\.4', "analysis/jets.py")
    require(jets, r'DeltaR\(\{p4_branch\}\[i\], mu2_p4\) > 0\.4', "analysis/jets.py")

    gen_filter = read(repo / "common/gen_vbf_filter.py")
    function = gen_filter[gen_filter.index("def ApplyGenVBFFilter"):]
    if ".Filter(" in function:
        raise RuntimeError("ApplyGenVBFFilter now applies an event filter: update its table row")

    require(skim, r'ApplyJetVetoMap\(df, config, muon_pt_default_suffix, False,', "analysis/skim.py")
    # muons_selection is evaluated in the skim only to store flag columns.
    function = muons[muons.index("def DefineMuonSelection"):]
    function = function[:function.index("\ndef ") if "\ndef " in function else len(function)]
    if ".Filter(" in function:
        raise RuntimeError("DefineMuonSelection now filters events: update the scope note")


def esc(value: Any) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def tf(value: Any) -> str:
    return "yes" if bool(value) else "no"


def make_rows(configs: dict[str, dict[str, Any]]) -> list[tuple[str, str, str, str, str, str]]:
    first = configs[ERAS[0]]["selections"]
    def common_setting(key: str, default: Any) -> str:
        values = {group["selections"].get(key, default) for group in configs.values()}
        return str(next(iter(values))) if len(values) == 1 else f"era-specific {key} (see era table)"

    mu_pt = common_setting("muon_pt_min", 15.0)
    mass_lo = common_setting("dimuon_mass_min", 50)
    mass_hi = common_setting("dimuon_mass_max", 200)
    jet_pt = common_setting("jet_pt_min", 20.0)
    jet_eta = common_setting("jet_eta_max", 4.7)
    met_flags = first_config_flags = configs[ERAS[0]]["main"].get("MET_flags", [])
    if not all(c["main"].get("MET_flags", []) == first_config_flags for c in configs.values()):
        met_text = "Era-specific logical AND: see config/Run3_*/maincfg.yaml"
    else:
        met_text = "Logical AND of " + ", ".join(flag.removeprefix("Flag_") for flag in met_flags)

    return [
        ("1", "Orthogonal luminosity split", "Event filter", "OrthogonalEraTag == 0/1/2 (seed 12345; luminosity fractions 110/111/26)", "MC 2024--2026", "Before all physics corrections"),
        ("2", "Generator VBF phase space", "Definition only", "Compute GenVBFFilter for selected inclusive/VBF-filtered DY samples; current code does not call df.Filter", "MC DY 2024--2026", "Before corrections; no rejection"),
        ("3", "Base event weights", "Weight definition", "Luminosity, generator sign, pileup and cross-section weights; configured QCD-scale sums", "MC", "Pileup is applied here; no rejection"),
        ("4", "Golden JSON", "Event filter", "Accept certified run:luminosityBlock pairs from maincfg.yaml", "Data", "Before object corrections"),
        ("5", "Muon ScaRe", "Object correction", "Momentum scale for data/MC; MC resolution smearing; NanoAOD and beam-spot-constrained momenta", "Data and MC", "Before every muon selection"),
        ("6", "Muon FSR recovery", "Object correction", "FSR photon recovery implemented in corrections/muon_fsr.py", "Data and MC", "After ScaRe; nominal selection uses Muon_pt_FSR_corr"),
        ("7", "JEC/JER", "Object correction", "Reapply JEC from the raw pT, Jet_pt * (1 - Jet_rawFactor); apply JER to MC; produce JER/JES Total shifts when variations are enabled", "Data and MC", "Before MET and jet selections; Jet_pt_nocorr keeps the NanoAOD pT, Jet_rawFactor is not stored"),
        ("8", "MET quality flags", "Event filter", met_text, "Data and MC", "NanoAOD flags; special 2022/2023 data bad-MET redefinition"),
        ("9", "Single-muon trigger matching", "Event filter", "HLT_IsoMu24; uncorrected muon pT (Muon_pt_raw_noCorr) > 26 GeV; trigger-object ID 13, filter bit 8, DeltaR(uncorrected muon, object) < 0.4", "Data and MC", "Uncorrected pT and direction (no ScaRe, no FSR), also for the muon shifts; OR over nominal/shifts"),
        ("10", "Good muons", "Object selection", f"uncorrected pT (Muon_pt_raw_noCorr) > {mu_pt} GeV; NanoAOD |eta| < 2.4; medium ID; pfIsoId >= 2", "Data and MC", "Cut on the uncorrected pT; the two muons are ordered by the ScaRe+FSR pT"),
        ("11", "Exactly two muons", "Event filter", "Exactly two good muons; nominal and every muon shift must pass when variations are enabled", "Data and MC", "After ScaRe, FSR and trigger matching"),
        ("12", "Dimuon mass window", "Event filter", f"{mass_lo} < m_mumu < {mass_hi} GeV; nominal and every muon shift must pass", "Data and MC", "Corrected mass: selected ScaRe+FSR muons"),
        ("13", "Electron veto", "Event filter", "No electron with pT > 20 GeV, |eta| < 2.5 and Electron_mvaIso_WP90", "Data and MC", "Raw NanoAOD electron kinematics"),
        ("14", "Muon ID/ISO/trigger SF", "Weight definition", "Configured correctionlib SF sources for both selected muons; central and optional up/down", "MC", "SFs evaluated at the uncorrected pT (pt_raw_noCorr) and the ScaRe+FSR eta; the tight/medium choice below 30 GeV uses the ScaRe+FSR pT"),
        ("15", "Jet ID and b-tag flags", "Object definition", "Tight jet ID; b-tag flags at |eta| < 2.5; b-tag event-veto flag is stored but not filtered", "Data and MC", "JEC/JER jets and era-specific WPs"),
        ("16", "Jet veto map", "Object definition", "Evaluate map; apply_filter=False in analysis/skim.py", "Data and MC", "Corrected jets; no direct event filter"),
        ("17", "Selected jets", "Object selection", f"corrected pT > {jet_pt} GeV; |eta| < {jet_eta}; tight ID; outside veto map; DeltaR(j,mu1/2) > 0.4", "Data and MC", "JEC/JER jets; repeated for JER/JES shifts"),
        ("18", "VBF candidates", "Definition only", "Build pair outside the horn; store HasVBF and indices", "Data and MC", "No VBF event filter in skim"),
        ("19", "Snapshot", "Output", "Write surviving events and nominal/shifted columns", "Data and MC", "Final step; categories are downstream"),
    ]


def render(repo: Path) -> str:
    validate_workflow(repo)
    configs = {
        era: {
            "main": load_yaml(repo / "config" / era / "maincfg.yaml"),
            "selections": load_yaml(repo / "config" / era / "selections.yaml"),
            "triggers": load_yaml(repo / "config" / era / "triggers.yaml"),
        }
        for era in ERAS
    }
    rows = make_rows(configs)
    lines = [
        r"% Generated by sync/python/collect_skim_selections.py; do not edit manually.",
        r"\begin{longtable}{>{\centering\arraybackslash}p{0.035\textwidth} >{\raggedright\arraybackslash}p{0.115\textwidth} >{\raggedright\arraybackslash}p{0.095\textwidth} >{\raggedright\arraybackslash}p{0.36\textwidth} >{\raggedright\arraybackslash}p{0.105\textwidth} >{\raggedright\arraybackslash}p{0.225\textwidth}}",
        r"\caption{Selections and corrections applied during NanoAOD skimming, in execution order.}\label{tab:skim-selections} \\",
        r"\toprule",
        r"\textbf{Step} & \textbf{Operation} & \textbf{Kind} & \textbf{Condition / definition} & \textbf{Scope} & \textbf{Correction state} \\",
        r"\midrule", r"\endfirsthead",
        r"\multicolumn{6}{c}{\tablename\ \thetable\ -- continued} \\",
        r"\toprule",
        r"\textbf{Step} & \textbf{Operation} & \textbf{Kind} & \textbf{Condition / definition} & \textbf{Scope} & \textbf{Correction state} \\",
        r"\midrule", r"\endhead",
        r"\midrule \multicolumn{6}{r}{Continued on next page} \\", r"\endfoot",
        r"\bottomrule", r"\endlastfoot",
    ]
    lines.extend(" & ".join(esc(cell) for cell in row) + r" \\" for row in rows)
    lines += [r"\end{longtable}", "", r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.11\textwidth} >{\raggedright\arraybackslash}p{0.18\textwidth} >{\raggedright\arraybackslash}p{0.31\textwidth} >{\raggedright\arraybackslash}p{0.34\textwidth}}", r"\caption{Era-dependent skim settings read from the YAML configuration.}\label{tab:skim-selection-eras} \\", r"\toprule", r"\textbf{Era} & \textbf{Muon/jet settings} & \textbf{Horn veto expression} & \textbf{Other event-filter settings} \\", r"\midrule", r"\endfirsthead", r"\toprule", r"\textbf{Era} & \textbf{Muon/jet settings} & \textbf{Horn veto expression} & \textbf{Other event-filter settings} \\", r"\midrule", r"\endhead", r"\bottomrule", r"\endlastfoot"]
    for era, group in configs.items():
        sel, main, trg = group["selections"], group["main"], group["triggers"]
        settings = f"mu pT>{sel.get('muon_pt_min')} GeV, jet pT>{sel.get('jet_pt_min')} GeV; |eta_mu|<2.4, |eta_j|<{sel.get('jet_eta_max')}"
        trigger_paths = [p for value in trg.values() if isinstance(value, dict) for p in value.get("path", [])]
        special = f"trigger filter={tf(sel.get('apply_trg_filter', True))}; jet-veto event filter={tf(sel.get('apply_jetveto_filter', False))}; HLT={','.join(trigger_paths)}"
        if main.get("badMET_flag_runs"):
            special += "; special bad-MET runs=" + "--".join(map(str, main["badMET_flag_runs"]))
        lines.append(" & ".join(esc(x) for x in (era, settings, sel.get("jet_horn_veto_expr", "false"), special)) + r" \\")
    lines += [r"\end{longtable}", "", r"\paragraph{Important scope note.}", r"The \texttt{muons\_selection} expressions in \texttt{config/Run3\_*/selections.yaml} are evaluated by \texttt{analysis/skim.py} only to store flag columns (for example \texttt{mu1\_pt\_g26}, \texttt{muons\_OS}); they reject no event. The \texttt{categories} and \texttt{masses\_regions} expressions are not invoked by the skim: they are downstream histogram-level definitions, not NanoAOD skim event filters.", ""]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="output .tex path")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    output = args.output or repo / "sync/latex/skim_selections_table.tex"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render(repo), encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
