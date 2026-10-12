
#!/usr/bin/env python3
import ROOT


ROOT.gInterpreter.Declare(
    """
    using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
    using RVecLV = ROOT::VecOps::RVec<LorentzVectorM>;
    using RVecF = ROOT::VecOps::RVec<float>;
    using RVecI = ROOT::VecOps::RVec<int>;

    LorentzVectorM fsr_corrected_p4(float mu_pt,
                                        float mu_eta,
                                        float mu_phi,
                                        float mu_mass,
                                        int mu_fsrIdx,
                                        const RVecF& fsr_pt,
                                        const RVecF& fsr_eta,
                                        const RVecF& fsr_phi,
                                        const RVecF& fsr_dROverEt2,
                                        const RVecF& fsr_relIso03,
                                        const RVecF& fsr_electronIdx) {
            LorentzVectorM res{mu_pt, mu_eta, mu_phi, mu_mass};

            if (mu_fsrIdx == -1) {
                return res;
            }

            float deltaR_mu_fsr = ROOT::VecOps::DeltaR(mu_eta, fsr_eta[mu_fsrIdx], mu_phi, fsr_phi[mu_fsrIdx]);

            if (!((deltaR_mu_fsr > 0.0001) && (deltaR_mu_fsr < 0.5))) {
                return res;
            }

            if (fsr_electronIdx[mu_fsrIdx] != -1) {
                return res;
            }

            if (fsr_pt[mu_fsrIdx] / mu_pt > 0.4) {
                return res;
            }

            if (fsr_dROverEt2[mu_fsrIdx] > 0.012) {
                return res;
            }

            if (fsr_relIso03[mu_fsrIdx] / mu_pt > 1.8) {
                return res;
            }

            res += LorentzVectorM{fsr_pt[mu_fsrIdx], fsr_eta[mu_fsrIdx], fsr_phi[mu_fsrIdx], 0.0};

            return res;
        }

        RVecLV fsr_corrected_p4(const RVecF& mu_pt,
                                const RVecF& mu_eta,
                                const RVecF& mu_phi,
                                const RVecF& mu_mass,
                                const RVecI& mu_fsrIdx,
                                const RVecF& fsr_pt,
                                const RVecF& fsr_eta,
                                const RVecF& fsr_phi,
                                const RVecF& fsr_dROverEt2,
                                const RVecF& fsr_relIso03,
                                const RVecF& fsr_electronIdx) {
            RVecLV res_vec(mu_pt.size());
            for (size_t i = 0; i < mu_pt.size(); ++i) {
                res_vec[i] = fsr_corrected_p4(mu_pt[i],
                                              mu_eta[i],
                                              mu_phi[i],
                                              mu_mass[i],
                                              mu_fsrIdx[i],
                                              fsr_pt,
                                              fsr_eta,
                                              fsr_phi,
                                              fsr_dROverEt2,
                                              fsr_relIso03,
                                              fsr_electronIdx);
            }
            return res_vec;
        }
        """
)

def apply_muon_fsr(df, is_data, want_variations):
    branches_map = {
        # uncorrected
        'Muon_pt':['Muon_p4_nano_FSR','Muon_pt_nano_FSR'],
        'Muon_bsConstrainedPt':['Muon_p4_bsc_FSR','Muon_pt_bsc_FSR'],
        # only scale
        'Muon_pt_nano_scale':['Muon_p4_nano_scale_FSR','Muon_pt_nano_scale_FSR'],
        'Muon_pt_bsc_scale':['Muon_p4_bsc_scale_FSR','Muon_pt_bsc_scale_FSR'],
        # scale & resolution
        'Muon_pt_nano_corr':['Muon_p4_nano_corr_FSR','Muon_pt_nano_corr_FSR'],
        'Muon_pt_bsc_corr':['Muon_p4_bsc_corr_FSR','Muon_pt_bsc_corr_FSR'],
    }
    branches_map_varied = {
        # scale & respolution varied
        'Muon_pt_nano_scale_up':['Muon_p4_nano_FSR_scale_up','Muon_pt_nano_FSR_scale_up'],
        'Muon_pt_nano_scale_down':['Muon_p4_nano_FSR_scale_down','Muon_pt_nano_FSR_scale_down'],
        'Muon_pt_nano_res_up':['Muon_p4_nano_FSR_res_up','Muon_pt_nano_FSR_res_up'],
        'Muon_pt_nano_res_down':['Muon_p4_nano_FSR_res_down','Muon_pt_nano_FSR_res_down'],
        'Muon_pt_bsc_scale_up':['Muon_p4_bsc_FSR_scale_up','Muon_pt_bsc_FSR_scale_up'],
        'Muon_pt_bsc_scale_down':['Muon_p4_bsc_FSR_scale_down','Muon_pt_bsc_FSR_scale_down'],
        'Muon_pt_bsc_res_up':['Muon_p4_bsc_FSR_res_up','Muon_pt_bsc_FSR_res_up'],
        'Muon_pt_bsc_res_down':['Muon_p4_bsc_FSR_res_down','Muon_pt_bsc_FSR_res_down'],
    }
    if not is_data and want_variations:
        branches_map.update(branches_map_varied)

    for muon_pt_branch,muon_pt_fsr_branch in branches_map.items():
        # if muon_pt_branch not in df.GetColumnNames(): continue
        df = df.Define(muon_pt_fsr_branch[0],f"fsr_corrected_p4({muon_pt_branch},Muon_eta, Muon_phi, Muon_mass, Muon_fsrPhotonIdx, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi, FsrPhoton_dROverEt2, FsrPhoton_relIso03, FsrPhoton_electronIdx)")
        df = df.Define(muon_pt_fsr_branch[1],f""" RVecF muon_pt_with_FSR({muon_pt_fsr_branch[0]}.size(), 0.);
                       for (size_t i = 0 ; i < {muon_pt_fsr_branch[0]}.size(); i++){{
                            muon_pt_with_FSR[i] = {muon_pt_fsr_branch[0]}[i].Pt();
                       }}
                       return muon_pt_with_FSR;
                       """)

    return df
