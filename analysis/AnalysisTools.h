#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecS = ROOT::VecOps::RVec<size_t>;
using RVecUC = ROOT::VecOps::RVec<UChar_t>;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecB = ROOT::VecOps::RVec<bool>;
using RVecVecI = ROOT::VecOps::RVec<RVecI>;
using RVecLV = ROOT::VecOps::RVec<LorentzVectorM>;
using RVecSetInt = ROOT::VecOps::RVec<std::set<int>>;

enum class GenLeptonMatch : int { Electron = 1, Muon = 2, TauElectron = 3, TauMuon = 4, Tau = 5, NoMatch = 6 };

enum class Period : int {
    Run2_2016_HIPM = 1,
    Run2_2016 = 2,
    Run2_2017 = 3,
    Run2_2018 = 4,
    Run3_2022 = 5,
    Run3_2022EE = 6,
    Run3_2023 = 7,
    Run3_2023BPix = 8,
    Run3_2024 = 9,
    Run3_2025 = 10,
    Run3_2026 = 11
};

RVecS CreateIndexes(size_t vecSize) {
    RVecS i(vecSize);
    std::iota(i.begin(), i.end(), 0);
    return i;
}

template <typename V>
RVecI ReorderObjects(const V &varToOrder, const RVecI &indices, size_t nMax = std::numeric_limits<size_t>::max()) {
    RVecI ordered_indices = indices;
    std::sort(ordered_indices.begin(), ordered_indices.end(), [&](int a, int b) {
        return varToOrder.at(a) > varToOrder.at(b);
    });
    const size_t n = std::min(ordered_indices.size(), nMax);
    ordered_indices.resize(n);
    return ordered_indices;
}


inline LorentzVectorM GetP4(const RVecF &pt, const RVecF &eta, const RVecF &phi, const RVecF &mass, int idx) {
    return LorentzVectorM(pt[idx], eta[idx], phi[idx], mass[idx]);
}

inline LorentzVectorM GetP4(const RVecF &pt, const RVecF &eta, const RVecF &phi, double mass, int idx) {
    return LorentzVectorM(pt[idx], eta[idx], phi[idx], mass);
}

RVecLV GetP4(const RVecF &pt, const RVecF &eta, const RVecF &phi, const RVecF &mass) {
    RVecLV p4;
    p4.reserve(pt.size());
    for (size_t idx = 0; idx < pt.size(); idx++)
        p4.emplace_back(pt[idx], eta[idx], phi[idx], mass[idx]);
    return p4;
}

RVecLV GetP4(const RVecF &pt, const RVecF &eta, const RVecF &phi, const RVecF &mass, const RVecS &indices) {
    RVecLV p4;
    p4.reserve(indices.size());
    for (auto &idx : indices)
        p4.emplace_back(pt[idx], eta[idx], phi[idx], mass[idx]);
    return p4;
}

RVecB RemoveOverlaps(const RVecLV &obj_p4,
                     const RVecB &pre_sel,
                     const std::vector<RVecLV> &other_objects,
                     size_t min_number_of_non_overlaps,
                     double min_deltaR) {
    RVecB result(pre_sel);
    const double min_deltaR2 = std::pow(min_deltaR, 2);

    const auto hasMinNumberOfNonOverlaps = [&](const LorentzVectorM &p4) {
        size_t cnt = 0;
        for (const auto &other_obj_col : other_objects) {
            for (const auto &other_obj_p4 : other_obj_col) {
                const double dR2 = ROOT::Math::VectorUtil::DeltaR2(p4, other_obj_p4);
                if (dR2 > min_deltaR2) {
                    ++cnt;
                    if (cnt >= min_number_of_non_overlaps)
                        return true;
                }
            }
        }
        return false;
    };

    for (size_t obj_idx = 0; obj_idx < obj_p4.size(); ++obj_idx) {
        result[obj_idx] = pre_sel[obj_idx] && hasMinNumberOfNonOverlaps(obj_p4.at(obj_idx));
    }
    return result;
}

RVecB RemoveOverlaps(const RVecLV &obj_p4, const RVecB &pre_sel, const RVecLV &other_objects, double min_deltaR) {
    RVecB result(pre_sel);
    const double min_deltaR2 = std::pow(min_deltaR, 2);

    const auto hasOverlaps = [&](const LorentzVectorM &p4) {
        for (const auto &other_obj_p4 : other_objects) {
            const double dR2 = ROOT::Math::VectorUtil::DeltaR2(p4, other_obj_p4);
            if (dR2 <= min_deltaR2)
                return true;
        }
        return false;
    };

    for (size_t obj_idx = 0; obj_idx < obj_p4.size(); ++obj_idx) {
        result[obj_idx] = pre_sel[obj_idx] && !hasOverlaps(obj_p4.at(obj_idx));
    }
    return result;
}

template <typename LVec, typename LVecCollection>
double MinDeltaR(const LVec &obj_p4, const LVecCollection &other_objects) {
    double min_dR = std::numeric_limits<double>::infinity();
    for (const auto &other_obj_p4 : other_objects) {
        const double dR = ROOT::Math::VectorUtil::DeltaR(obj_p4, other_obj_p4);
        if (dR < min_dR) {
            min_dR = dR;
        }
    }
    return min_dR;
}

int FindMatching(const LorentzVectorM &target_p4, const RVecLV &ref_p4, const float deltaR_thr) {
    double deltaR_min = deltaR_thr;
    int current_idx = -1;
    for (int refIdx = 0; refIdx < ref_p4.size(); refIdx++) {
        auto dR_targetRef = ROOT::Math::VectorUtil::DeltaR(target_p4, ref_p4.at(refIdx));
        if (dR_targetRef < deltaR_min) {
            deltaR_min = dR_targetRef;
            current_idx = refIdx;
        }
    }
    return current_idx;
}

RVecI FindMatching(const RVecLV &target_p4, const RVecLV &ref_p4, const float deltaR_thr) {
    RVecI targetIndices(target_p4.size(), -1);
    for (int targetIdx = 0; targetIdx < target_p4.size(); targetIdx++) {
        int refIdxFound = FindMatching(target_p4[targetIdx], ref_p4, deltaR_thr);
        targetIndices[targetIdx] = refIdxFound;
    }
    return targetIndices;
}

int FindMatching(const bool pre_sel_target,
                 const RVecB &pre_sel_ref,
                 const LorentzVectorM &target_p4,
                 const RVecLV &ref_p4,
                 const float dR_thr) {
    // RVecI matched(1,-1); // Only one target, so size is 1 and initialized with
    // false
    int current_idx = -1;
    float deltaR_min = dR_thr;
    if (pre_sel_target) {
        for (size_t ref_idx = 0; ref_idx < pre_sel_ref.size(); ref_idx++) {
            if (pre_sel_ref[ref_idx] == 0)
                continue;
            auto dR_targetRef = ROOT::Math::VectorUtil::DeltaR(target_p4, ref_p4[ref_idx]);
            if (dR_targetRef < deltaR_min) {
                current_idx = ref_idx;
                deltaR_min = dR_targetRef;
            }
        }
    }
    return current_idx;
}

RVecI FindMatching(const RVecB &pre_sel_target,
                   const RVecB &pre_sel_ref,
                   const RVecLV &target_p4,
                   const RVecLV &ref_p4,
                   const float deltaR_thr) {
    RVecI targetIndices(target_p4.size(), -1);
    for (int targetIdx = 0; targetIdx < target_p4.size(); targetIdx++) {
        targetIndices[targetIdx] =
            FindMatching(pre_sel_target[targetIdx], pre_sel_ref, target_p4[targetIdx], ref_p4, deltaR_thr);
    }
    return targetIndices;
}

RVecSetInt FindMatchingSet(const RVecB &pre_sel_target,
                           const RVecB &pre_sel_ref,
                           const RVecLV &target_p4,
                           const RVecLV &ref_p4,
                           const float dR_thr) {
    RVecSetInt findMatching(pre_sel_target.size());
    for (size_t ref_idx = 0; ref_idx < pre_sel_ref.size(); ref_idx++) {
        if (pre_sel_ref[ref_idx] == 0)
            continue;
        for (size_t target_idx = 0; target_idx < pre_sel_target.size(); target_idx++) {
            if (pre_sel_target[target_idx] == 0)
                continue;
            auto dR_current = ROOT::Math::VectorUtil::DeltaR(target_p4[target_idx], ref_p4[ref_idx]);
            if (dR_current < dR_thr) {
                findMatching[target_idx].insert(ref_idx);
            }
        }
    }
    return findMatching;
}

namespace ROOT {
    namespace VecOps {
        template <typename TIn, typename TOut>
        RVec<TOut> TakeAndCast(const RVec<TIn> &v,
                               const RVec<typename RVec<TIn>::size_type> &i,
                               const TOut default_val) {
            RVec<TOut> result(i.size());
            for (typename RVec<TIn>::size_type pos = 0; pos < i.size(); ++pos) {
                const auto idx = i[pos];
                result[pos] = idx >= 0 && idx < v.size() ? static_cast<TOut>(v[idx]) : default_val;
            }
            return result;
        }
    }  // namespace VecOps
}  // namespace ROOT

namespace v_ops {
    template <typename VecIn, typename OutT, auto MemPtr>
    ROOT::VecOps::RVec<OutT> extract(const VecIn &p4) {
        return ROOT::VecOps::Map(p4, [](const auto &v) -> OutT { return static_cast<OutT>((v.*MemPtr)()); });
    }

#define DEFINE_EXTRACTOR(fn_name, method_name)                      \
    template <typename LV, typename OutT = float>                   \
    ROOT::VecOps::RVec<OutT> fn_name(const LV &p4) {                \
        return extract<LV, OutT, &LV::value_type::method_name>(p4); \
    }

    DEFINE_EXTRACTOR(pt, pt)
    DEFINE_EXTRACTOR(eta, eta)
    DEFINE_EXTRACTOR(phi, phi)
    DEFINE_EXTRACTOR(mass, mass)
    DEFINE_EXTRACTOR(energy, energy)
    DEFINE_EXTRACTOR(Et, Et)
    DEFINE_EXTRACTOR(rapidity, Rapidity)

#undef DEFINE_EXTRACTOR
}  // namespace v_ops

namespace eventId {

    ULong64_t encodeFullEventId(ULong64_t sample_name_crc, ULong64_t infile_crc, ULong64_t rdfentry) {
        if (sample_name_crc >> 16)
            throw std::runtime_error("sample_name_crc overflows 16 bits");
        if (infile_crc >> 16)
            throw std::runtime_error("infile_crc overflows 16 bits");
        if (rdfentry >> 32)
            throw std::runtime_error("rdfentry overflows 32 bits");

        return (sample_name_crc << 48) | (infile_crc << 32) | rdfentry;
    }

    std::tuple<ULong64_t, ULong64_t, ULong64_t> decodeFullEventId(ULong64_t fullEventId) {
        ULong64_t sample_name_crc = (fullEventId >> 48) & 0xFFFF;
        ULong64_t infile_crc = (fullEventId >> 32) & 0xFFFF;
        ULong64_t rdfentry = fullEventId & 0xFFFFFFFF;

        return std::make_tuple(sample_name_crc, infile_crc, rdfentry);
    }

}  // namespace eventId




struct VBFJets {
    static constexpr size_t n_legs = 2;
    static constexpr float m_inv_th = 400.;
    static constexpr float eta_th = 2.5;
    static constexpr float leading_pt_th = 35;
    static constexpr float subleading_pt_th = 25;
    std::array<int, n_legs> leg_index;
    std::array<LorentzVectorM, n_legs> leg_p4;  // p4 of quark from H->bb
    float m_inv;
    float eta_separation;
    bool isVBF;
    RVecLV legs_p4;
};

VBFJets FindVBFJets(const RVecLV& Jet_p4, const RVecB& pre_sel) {
    VBFJets VBF_jets_collection;
    // initialize with false value
    VBF_jets_collection.isVBF = false;
    // initialize with negative indices
    for (int j = 0; j < VBF_jets_collection.n_legs; j++) {
        VBF_jets_collection.leg_index[j] = -1;
    }
    // if (Jet_p4.size() < VBF_jets_collection.n_legs)
    //     return VBF_jets_collection;

    float inv_mass_th = VBF_jets_collection.m_inv_th;
    float eta_th = VBF_jets_collection.eta_th;
    float leading_pt_th = VBF_jets_collection.leading_pt_th;
    float subleading_pt_th = VBF_jets_collection.subleading_pt_th;

    for (size_t j1_idx = 0; j1_idx < Jet_p4.size(); j1_idx++) {
        if (pre_sel[j1_idx]==0)
            continue;
        for (size_t j2_idx = j1_idx +1 ; j2_idx < Jet_p4.size(); j2_idx++) {  // j > i to avoid duplicates
            if (pre_sel[j2_idx]==0)
                continue;
            // comparison with Run2: same selection. Since there could be MORE than one jet pair, the one with highest invariant mass is selected, keeping the DeltaEta threshold fixed to 2.5
            float inv_mass = (Jet_p4.at(j1_idx) + Jet_p4.at(j2_idx)).M();
            float eta = Jet_p4.at(j1_idx).Eta() - Jet_p4.at(j2_idx).Eta();
            float leadingjet_pt = Jet_p4.at(j1_idx).Pt();
            float subleadingjet_pt = Jet_p4.at(j2_idx).Pt();
            if (inv_mass >= inv_mass_th && std::abs(eta) >= eta_th && leadingjet_pt>=leading_pt_th && subleadingjet_pt>=subleading_pt_th) {
                inv_mass_th = inv_mass;
                // eta_th = eta;
                VBF_jets_collection.leg_index[0] = j1_idx;
                VBF_jets_collection.leg_index[1] = j2_idx;
                VBF_jets_collection.leg_p4[0] = Jet_p4.at(j1_idx);
                VBF_jets_collection.leg_p4[1] = Jet_p4.at(j2_idx);
                VBF_jets_collection.m_inv = inv_mass;
                VBF_jets_collection.eta_separation = eta;
                VBF_jets_collection.isVBF = true;
                VBF_jets_collection.legs_p4 = {Jet_p4.at(j1_idx), Jet_p4.at(j2_idx)};
            }
        }
    }
    return VBF_jets_collection;
}

std::pair<double, double> ComputeCosThetaPhiCS(const LorentzVectorM& mu1_p4,
                                               const LorentzVectorM& mu2_p4,
                                               double Ebeam) {
    // muons p4 in XYZ coordinates
    LorentzVectorXYZ mu1_p4_XYZ = LorentzVectorXYZ{mu1_p4.Px(), mu1_p4.Py(), mu1_p4.Pz(), mu1_p4.E()};
    LorentzVectorXYZ mu2_p4_XYZ = LorentzVectorXYZ{mu2_p4.Px(), mu2_p4.Py(), mu2_p4.Pz(), mu2_p4.E()};
    // dilepton boosted p4 in XYZ coordinates
    LorentzVectorXYZ dilepton = mu1_p4_XYZ + mu2_p4_XYZ;
    // boost vector in XYZ coordinates
    ROOT::Math::XYZVector boost = -dilepton.BoostToCM();
    // boost muons in XYZ coordinates
    LorentzVectorXYZ mu1_boosted_XYZ = ROOT::Math::VectorUtil::boost(mu1_p4_XYZ, boost);
    LorentzVectorXYZ mu2_boosted_XYZ = ROOT::Math::VectorUtil::boost(mu2_p4_XYZ, boost);
    // proton p4 in XYZ coordinates
    LorentzVectorXYZ pA(0, 0, Ebeam, Ebeam);
    LorentzVectorXYZ pB(0, 0, -Ebeam, Ebeam);
    // boost proton p4 in XYZ coordinates
    LorentzVectorXYZ pA_boosted = ROOT::Math::VectorUtil::boost(pA, boost);
    LorentzVectorXYZ pB_boosted = ROOT::Math::VectorUtil::boost(pB, boost);
    // Collins-Soper axes
    // z axis is the bisector of the angle between the two boosted muons
    // y axis is the normal to the plane defined by the two boosted muons and the two boosted protons
    const ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> z_cs =
        (pA_boosted.Vect().Unit() - pB_boosted.Vect().Unit()).Unit();
    const ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> y_cs =
        (pA_boosted.Vect().Unit().Cross(pB_boosted.Vect().Unit())).Unit();
    const ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> x_cs = (y_cs.Cross(z_cs)).Unit();
    // cos(theta_CS)
    double cos_theta_cs = mu1_boosted_XYZ.Vect().Unit().Dot(z_cs);
    // phi_CS
    double phi_cs = atan2(mu1_boosted_XYZ.Vect().Dot(y_cs), mu1_boosted_XYZ.Vect().Dot(x_cs));

    return std::make_pair(cos_theta_cs, phi_cs);
}



float pT_sum(const RVecLV& all_p4s) {
    float pT_x_sum = 0.;
    float pT_y_sum = 0.;
    for (int idx = 0; idx < all_p4s.size(); idx++) {
        pT_x_sum += all_p4s[idx].Px();
        pT_y_sum += all_p4s[idx].Py();
    }
    return sqrt(pT_x_sum * pT_x_sum + pT_y_sum * pT_y_sum);
}

float pT_diff(const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
    float delta_px = p4_1.Px() - p4_2.Px();
    float delta_py = p4_1.Py() - p4_2.Py();
    return sqrt(delta_px * delta_px + delta_py * delta_py);
}
