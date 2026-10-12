#!/usr/bin/env python3
"""Generate the compact, analyst-facing Run 3 selection table."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml


ERAS = ("2022", "2022EE", "2023", "2023BPix", "2024", "2025", "2026")


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return yaml.safe_load(stream) or {}


def esc(value: Any) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def common(configs: list[dict[str, Any]], key: str) -> Any:
    values = {config.get(key) for config in configs}
    if len(values) != 1:
        raise ValueError(f"selection setting {key!r} differs across eras: {values}")
    return values.pop()


def render(repo: Path) -> str:
    selections = [load(repo / "config" / f"Run3_{era}" / "selections.yaml") for era in ERAS]
    main_configs = [load(repo / "config" / f"Run3_{era}" / "maincfg.yaml") for era in ERAS]
    pt_mu = common(selections, "muon_pt_min")
    mass_min = common(selections, "dimuon_mass_min")
    mass_max = common(selections, "dimuon_mass_max")
    pt_jet = common(selections, "jet_pt_min")
    eta_jet = common(selections, "jet_eta_max")
    taggers = [config.get("bTagAlgo", "PNet") for config in main_configs]
    if taggers[:4] != ["PNet"] * 4 or taggers[4:] != ["UParTAK4"] * 3:
        raise ValueError(f"unexpected era-dependent b-tag algorithms: {taggers}")

    rows = [
        (1, "Golden JSON", "Data only: require the run and luminosity section to be included in the era-specific Golden JSON certification mask.", r"$\ast$"),
        (2, "Orthogonal MC era split", r"MC only in 2024--2026: assign each event deterministically to exactly one era using run, luminosity block, event and dataset identifiers (seed 12345). The fractions follow the target luminosities 110:111:26 for 2024:2025:2026, preventing event reuse across eras.", r"$\ast$"),
        (3, "Generator VBF phase space", r"MC only in 2024--2026, for DY $105<m_{\mu\mu}<160$ samples. Remove hard-process leptons within $\Delta R<0.3$ of GenJets and define GenVBFFilter from the two leading remaining GenJets: pass when $m_{jj}>300$ GeV. The inclusive samples DYto2Mu\_MLL105To160 and DYto2Mu\_MLL\_105to160\_amcatnloFXFX use GenVBFFilter$=0$; their FlashSim aliases follow the same rule. DYto2Mu\_MLL105To160\_VBFFiltered and DYto2Mu\_MLL\_105to160\_amcatnloFXFX\_Fil\_VBF use GenVBFFilter$=1$. The flag is stored at skim level and the complementary cuts are applied during histogram routing.", r"$\ast$ / $\dagger$"),
        (4, "Event quality", "Logical AND of the recommended era-dependent NanoAOD MET filters.", r"$\ast$"),
        (5, "Muon momentum", r"For each muon use the beam-spot-constrained (BSC) $p_T$ when its BSC fit $\chi^2<30$; otherwise use the standard NanoAOD $p_T$. ScaRe scale corrections are applied to data and MC, resolution smearing to MC, and FSR recovery is included before the invariant-mass construction and the analysis $p_T$ cuts. The skim muon $p_T$ cut, the trigger matching and the muon scale factors use the uncorrected $p_T$ (no ScaRe, no FSR).", r"$\ast$"),
        (6, "Single-muon trigger", r"Require HLT\_IsoMu24 online and at least one selected offline muon with uncorrected $p_T>26$ GeV matched to a trigger object with ID 13, filter bit 8 and $\Delta R<0.4$ (uncorrected muon direction).", r"$\ast$ / $\dagger$"),
        (7, "Muon skim selection", rf"Exactly two muons with uncorrected $p_T>{pt_mu}$ GeV (Nano/BSC choice above, before ScaRe and FSR), $|\eta|<2.4$, medium ID and pfIsoId $\geq2$; they are ordered by the ScaRe- and FSR-corrected $p_T$.", r"$\ast$"),
        (8, "Dimuon skim window", rf"${mass_min}<m_{{\mu\mu}}<{mass_max}$ GeV, evaluated with the selected corrected muon momenta.", r"$\ast$"),
        (9, "Electron veto", r"No electron with $p_T>20$ GeV, $|\eta|<2.5$ and MVA Iso WP90.", r"$\ast$"),
        (10, "Muon analysis selection", r"Opposite-sign muons; corrected $p_T>26$ GeV for the leading and $>20$ GeV for the subleading muon, medium ID and pfIsoId $\geq2$. The online HLT\_IsoMu24 decision and the requirement that at least one of the two offline muons is trigger-matched are retained.", r"$\dagger$"),
        (11, "Low-pT muon quality", r"A selected muon with corrected $p_T<30$ GeV must additionally satisfy tight ID and pfIsoId $\geq4$ when the corresponding category is used.", r"$\dagger$"),
        (12, "Jet selection", rf"Jets with $p_T>{pt_jet}$ GeV after jet-energy scale and resolution corrections, $|\eta|<{eta_jet}$, tight jet ID, outside the era-dependent detector horn/veto region and $\Delta R(j,\mu)>0.4$ from both selected muons.", r"$\ast$"),
        (13, "b-jet veto", r"Within $|\eta|<2.5$, require zero medium-tagged jets and at most one loose-tagged jet. Working points are era dependent; the tagger is ParticleNet (PNet) for 2022--2023BPix and UParTAK4 for 2024--2026.", r"$\dagger$"),
        (14, "Mass regions", r"Inclusive: $60<m_{\mu\mu}<150$ GeV; Z sideband: $70$--$110$ GeV; signal fit: $115$--$135$ GeV; H sidebands: $110$--$115$ or $135$--$150$ GeV.", r"$\dagger$"),
        (15, "Baseline category", r"Full event selection (Golden JSON for data or MC-specific filters, MET quality and electron veto) plus the two-muon analysis selection, HLT\_IsoMu24, at least one offline muon matched to the fired online trigger, and the b-jet veto.", r"$\dagger$"),
        (16, "VBF topology and category", r"At skim level, among selected jets outside the horn, define the pair with the largest $m_{jj}$ satisfying $m_{jj}\geq400$ GeV and $|\Delta\eta_{jj}|\geq2.5$. At histogram level, the VBF category requires baseline and leading/subleading VBF-jet $p_T\geq35/25$ GeV.", r"$\ast_{\mathrm{def.}}$ / $\dagger$"),
        (17, "ggF categories", r"Baseline events failing the VBF definition, split into 0-jet, 1-jet and $\geq2$-jet categories when requested.", r"$\dagger$"),
    ]

    lines = [
        r"% Automatically generated analyst-facing summary.",
        r"\begin{longtable}{>{\centering\arraybackslash}p{0.045\textwidth} >{\raggedright\arraybackslash}p{0.19\textwidth} >{\raggedright\arraybackslash}p{0.64\textwidth} >{\centering\arraybackslash}p{0.07\textwidth}}",
        r"\caption{Run 3 analysis selection in application order.}\label{tab:analysis-selection} \\",
        r"\toprule",
        r"\textbf{Order} & \textbf{Requirement} & \textbf{Selection} & \textbf{Stage} \\",
        r"\midrule", r"\endfirsthead",
        r"\multicolumn{4}{c}{\tablename\ \thetable\ -- continued} \\",
        r"\toprule",
        r"\textbf{Order} & \textbf{Requirement} & \textbf{Selection} & \textbf{Stage} \\",
        r"\midrule", r"\endhead", r"\bottomrule", r"\endlastfoot",
    ]
    for order, requirement, selection, stage in rows:
        lines.append(f"{order} & {esc(requirement)} & {selection} & {stage} " + r"\\")
    lines += [
        r"\end{longtable}", "",
        r"\noindent\textbf{Stage notation:} $\ast$ applied during NanoAOD skimming; "
        r"$\dagger$ applied when producing analysis histograms. Requirements marked with both symbols "
        r"are enforced at both stages; $\ast_{\mathrm{def.}}$ means that the object or flag is defined but not used to reject events at skim level. Object calibrations and systematic variations are propagated before "
        r"the corresponding selections.", "",
        r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.15\textwidth} >{\raggedright\arraybackslash}p{0.23\textwidth} >{\raggedright\arraybackslash}p{0.55\textwidth}}",
        r"\caption{Reco-jet/GenJet component split for DY, EWK, ggF and VBF MC samples. This classification is orthogonal to the analysis category.}\label{tab:jet-gen-components} \\",
        r"\toprule",
        r"\textbf{Reco-jet multiplicity} & \textbf{Component} & \textbf{Definition} \\",
        r"\midrule", r"\endfirsthead",
        r"\toprule",
        r"\textbf{Reco-jet multiplicity} & \textbf{Component} & \textbf{Definition} \\",
        r"\midrule", r"\endhead", r"\bottomrule", r"\endlastfoot",
        r"Any & Matching rule & A selected reconstructed jet is hard-scatter matched when its NanoAOD GenJet index is non-negative. A negative index means that no GenJet is matched and the reconstructed jet is classified as a pileup (PU) jet. \\",
        r"0 jets & 0J Hard & No selected reconstructed jet; the event forms the exclusive zero-jet component. \\",
        r"1 jet & 1J Hard / 1J PU & Classify the leading selected reconstructed jet as hard or PU using its GenJet match. \\",
        r"$\geq2$ jets & 2J Hard / 2J PU1 / 2J PU2 & Count unmatched jets among the first two selected reconstructed jets: respectively 0, 1 or 2 PU jets. \\",
        r"VBF pair & VBF Hard / VBF PU1 / VBF PU2 & Apply the same 0/1/2 unmatched-jet count specifically to the two jets forming the selected VBF pair. \\",
        r"\end{longtable}", "",
        r"The generic 0/1/2 matching flags are process independent. Component outputs are produced for the configured DY/EWK/ggF/VBF process families when the jet--GenJet splitting workflow is enabled.", "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    output = args.output or repo / "sync/latex/analysis_selections_table.tex"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render(repo), encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
