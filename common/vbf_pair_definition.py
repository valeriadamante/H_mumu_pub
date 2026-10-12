"""Alternative definitions of the VBF jet pair, applied to the raw node.

All three require the basic VBF cuts (analysis/AnalysisTools.h,
struct VBFJets): m(jj) >= 400 GeV, |Delta eta(jj)| >= 2.5, pT >= 35/25 GeV.

  maxmjj   among the pairs of preselected jets passing the cuts, the one with
           the largest m(jj).  It is FindVBFJets, i.e. the production: no
           redefinition, the skim columns are kept.
  hardest  among all pairs of preselected jets passing the cuts, the hardest
           in the pT of the single jets: first the hardest leading, for equal
           leading the hardest subleading (lexicographic order in the
           SelectedJet collection, sorted in pT).
  leading  the two highest-pT jets among all reconstructed jets (pT > 25 GeV,
           separated from the muons): the pair exists only if both pass the
           preselection (ID, |eta|, veto map, horn) and the cuts.  A leading jet
           that falls in the horn removes the event from VBF instead of being
           replaced by the next one: this is the only difference with hardest.

Only HasVBF and VBFJetIdx_1/2 are overwritten: m_jj, vbfjet*_*, the DNN
variables and the Hard/PU DY components are recomputed from there
(common/add_vars.py, common/jet_component_splitting.py).  Only for the nominal
jet variation.
"""

from __future__ import annotations

import ROOT

DEFINITIONS = ("maxmjj", "hardest", "leading")

_DECLARED = False

_CODE = r"""
namespace vbfdef {

using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecB = ROOT::VecOps::RVec<bool>;
using P4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;

constexpr float kMjj = 400.f, kDeta = 2.5f, kPt1 = 35.f, kPt2 = 25.f, kJetPtMin = 25.f;

// Does the SelectedJet pair (i, j) pass the basic cuts? i is the harder jet.
inline bool PassesCuts(const RVecF& pt, const RVecF& eta, const RVecF& phi,
                       const RVecF& mass, int i, int j) {
    const P4 a(pt[i], eta[i], phi[i], mass[i]), b(pt[j], eta[j], phi[j], mass[j]);
    return (a + b).M() >= kMjj && std::abs(eta[i] - eta[j]) >= kDeta
        && pt[i] >= kPt1 && pt[j] >= kPt2;
}

// Indices {leg1, leg2} in the SelectedJet collection, sorted in pT; {-1,-1} if none.
RVecI Hardest(const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass) {
    // SelectedJet is sorted in pT: the first pair (i, j) that passes is the hardest.
    RVecI out{-1, -1};
    for (int i = 0; i < (int)pt.size(); ++i)
        for (int j = i + 1; j < (int)pt.size(); ++j)
            if (PassesCuts(pt, eta, phi, mass, i, j)) { out[0] = i; out[1] = j; return out; }
    return out;
}

RVecI Leading(const RVecF& jet_pt, const RVecF& jet_eta, const RVecF& jet_phi,
              const RVecB& jet_good, const RVecI& selected_idx,
              float mu1_eta, float mu1_phi, float mu2_eta, float mu2_phi,
              const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass) {
    RVecI out{-1, -1};
    int lead[2] = {-1, -1};
    for (int k = 0; k < (int)jet_pt.size(); ++k) {
        if (jet_pt[k] <= kJetPtMin) continue;
        if (ROOT::VecOps::DeltaR(jet_eta[k], mu1_eta, jet_phi[k], mu1_phi) <= 0.4f) continue;
        if (ROOT::VecOps::DeltaR(jet_eta[k], mu2_eta, jet_phi[k], mu2_phi) <= 0.4f) continue;
        if (lead[0] < 0 || jet_pt[k] > jet_pt[lead[0]]) { lead[1] = lead[0]; lead[0] = k; }
        else if (lead[1] < 0 || jet_pt[k] > jet_pt[lead[1]]) lead[1] = k;
    }
    if (lead[1] < 0 || !jet_good[lead[0]] || !jet_good[lead[1]]) return out;
    int sel[2] = {-1, -1};
    for (int s = 0; s < (int)selected_idx.size(); ++s)
        for (int l = 0; l < 2; ++l)
            if (selected_idx[s] == lead[l]) sel[l] = s;
    if (sel[0] < 0 || sel[1] < 0) return out;
    if (sel[0] > sel[1]) std::swap(sel[0], sel[1]);
    if (!PassesCuts(pt, eta, phi, mass, sel[0], sel[1])) return out;
    out[0] = sel[0]; out[1] = sel[1];
    return out;
}

}
"""

_EXPRESSION = {
    "hardest": "vbfdef::Hardest(SelectedJet_pt, SelectedJet_eta, SelectedJet_phi, SelectedJet_mass)",
    "leading": ("vbfdef::Leading(Jet_pt, Jet_eta, Jet_phi, Jet_NoOverlapWithMuons, SelectedJet_idx, "
                "mu1_eta, mu1_phi, mu2_eta, mu2_phi, "
                "SelectedJet_pt, SelectedJet_eta, SelectedJet_phi, SelectedJet_mass)"),
}


def apply_vbf_pair_definition(rdf, definition):
    """Redefine HasVBF and VBFJetIdx_1/2 on the raw skim node."""
    global _DECLARED
    if definition in (None, "", "maxmjj"):
        return rdf
    if definition not in _EXPRESSION:
        raise ValueError(f"Unknown VBF pair definition {definition!r}; available {DEFINITIONS}")
    if not _DECLARED:
        ROOT.gInterpreter.Declare(_CODE)
        _DECLARED = True
    rdf = rdf.Define("VBFPairIdx__def", _EXPRESSION[definition])
    rdf = rdf.Redefine("VBFJetIdx_1", "VBFPairIdx__def[0]")
    rdf = rdf.Redefine("VBFJetIdx_2", "VBFPairIdx__def[1]")
    return rdf.Redefine("HasVBF", "VBFPairIdx__def[0] >= 0")
