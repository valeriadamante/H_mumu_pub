
"""Physics observables and histogram selections added to dataframes."""


def _selection_suffixes(syst_cfg=None, want_variations=False):
    suffixes = [("", "", "")]
    if not want_variations or not syst_cfg:
        return suffixes
    scales = syst_cfg.get("scales", ["up", "down"])
    for syst_name, syst_info in syst_cfg.get("systematics", {}).items():
        if syst_name == "Central":
            continue
        components = syst_info.get("components", ())
        for component in components or (None,):
            for scale in scales:
                output_name = component or syst_name
                jet_suffix = syst_info.get("jet_suffix", "").format(scale=scale)
                if component:
                    jet_suffix = f"_{component}{scale.capitalize()}"
                suffixes.append((
                    f"_{output_name}{scale.capitalize()}",
                    syst_info.get("muon_suffix", "").format(scale=scale),
                    jet_suffix,
                ))
    return suffixes


def GetSelectionSuffixForSystematic(syst_name, syst_info=None):
    if syst_name == "Central" or syst_info is None:
        return ""
    if not syst_info.get("muon_suffix", "") and not syst_info.get("jet_suffix", ""):
        return ""
    return f"_{syst_name}"


def DefineSelections(df, sel_config, syst_cfg=None, want_variations=False):
    defined_columns = {str(column) for column in df.GetColumnNames()}
    section_suffix_key = {
        "masses_regions": "tot",
        "muons_selection": "mu",
        "jets_selection": "jet",
        "categories": "tot",
    }
    for section, suffix_key in section_suffix_key.items():
        for selection_name, content in sel_config.get(section, {}).items():
            expression_template = content.get("expression", "") if isinstance(content, dict) else content
            if not expression_template:
                print(f"[WARNING] Empty selection expression for {selection_name}. Skipping.")
                continue
            for total_suffix, muon_suffix, jet_suffix in _selection_suffixes(syst_cfg, want_variations):
                suffix = {"tot": total_suffix, "mu": muon_suffix, "jet": jet_suffix}[suffix_key]
                column_name = f"{selection_name}{suffix}"
                expression = expression_template.format(
                    tot_suff=total_suffix, mu_suff=muon_suffix, jet_suff=jet_suffix
                )
                # The selections are always recomputed from the expression, even
                # when a branch with the same name exists in the skim. Trusting
                # the stored branch breaks mixed-schema datasets: RDF derives the
                # columns from the first files of the chain, so a branch present
                # only in some files looks available and reading then fails on
                # the others with "tree does not have a branch".
                # Happened on Run3_2024/Muon1_Run2024G, where muons_SS_presel_trg
                # exists in 12 files out of 204.
                if column_name in defined_columns:
                    df = df.Redefine(column_name, expression)
                else:
                    df = df.Define(column_name, expression)
                    defined_columns.add(column_name)
    return df


def SelectedJetObservablesDef(df, suffix=""):
    columns = {str(col) for col in df.GetColumnNames()}
    jet_names = {
        0: "leading",
        1: "subleading",
        2: "third",
        3: "fourth",
    }
    sel_jet_cols = ["SelectedJet_area","SelectedJet_btagDeepFlavQG","SelectedJet_btagPNetB","SelectedJet_btagPNetQvG","SelectedJet_btagPNetCvL","SelectedJet_btagUParTAK4QvG","SelectedJet_eta","SelectedJet_idx","SelectedJet_mass","SelectedJet_phi","SelectedJet_pt"]
    for jet_idx, jet_type in jet_names.items():
        for jet_obs in sel_jet_cols:
            source = f"{jet_obs}{suffix}"
            if source not in columns:
                continue
            jet_obs_suff = "_".join(jet_obs.split("_")[1:])
            df = df.Define(
                f"{jet_type}jet_{jet_obs_suff}{suffix}",
                f"(SelectedJet_idx{suffix}.size()>{jet_idx} && !SelectedJet_IsInsideHorn{suffix}[{jet_idx}]) ? {source}.at({jet_idx}): -1000.f;",
            )
        df = df.Define(
            f"{jet_type}jet_p4{suffix}",
            f"(SelectedJet_idx{suffix}.size()>{jet_idx} && !SelectedJet_IsInsideHorn{suffix}[{jet_idx}]) ? ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(SelectedJet_pt{suffix}.at({jet_idx}), SelectedJet_eta{suffix}.at({jet_idx}),SelectedJet_phi{suffix}.at({jet_idx}), SelectedJet_mass{suffix}.at({jet_idx})) : ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(-10000.,-10000.,-10000.,-10000.);",
        )
    df = df.Define(
        f"delta_eta_jj_ls{suffix}",
        f"SelectedJet_idx{suffix}.size() >= 2 && std::abs(leadingjet_eta{suffix}) < 10. && std::abs(subleadingjet_eta{suffix}) < 10. ? std::abs(leadingjet_eta{suffix} - subleadingjet_eta{suffix}) : -1000.f"
    )
    df = df.Define(
        f"m_jj_ls{suffix}",
        f"leadingjet_p4{suffix}.M()>=0 && subleadingjet_p4{suffix}.M() >=0 ? (leadingjet_p4{suffix}+subleadingjet_p4{suffix}).M(): -1000.f",
    )

    return df


def VBFJetObservablesDef(df, suffix=""):
    columns = {str(col) for col in df.GetColumnNames()}
    # sel_jet_cols = [c for c in columns if c.startswith("SelectedJet")]
    sel_jet_cols = ["SelectedJet_area","SelectedJet_btagDeepFlavQG","SelectedJet_btagPNetB","SelectedJet_btagPNetQvG","SelectedJet_btagPNetCvNotB","SelectedJet_btagUParTAK4CvB","SelectedJet_btagUParTAK4QvG","SelectedJet_btagUParTAK4CvNotB","SelectedJet_btagUParTAK4B","SelectedJet_btagUParTAK4CvL","SelectedJet_btagUParTAK4UDG","SelectedJet_eta","SelectedJet_idx","SelectedJet_mass","SelectedJet_phi","SelectedJet_pt"]
    for vbfj_idx in [1,2]:
        df = df.Define(
            f"vbfjet{vbfj_idx}_p4{suffix}",
            f"(HasVBF{suffix}) ? ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(SelectedJet_pt{suffix}.at(VBFJetIdx_{vbfj_idx}{suffix}), SelectedJet_eta{suffix}.at(VBFJetIdx_{vbfj_idx}{suffix}),SelectedJet_phi{suffix}.at(VBFJetIdx_{vbfj_idx}{suffix}), SelectedJet_mass{suffix}.at(VBFJetIdx_{vbfj_idx}{suffix})) : ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(-10000.,-10000.,-10000.,-10000.);",
        )
        for jet_obs in sel_jet_cols:
            source = f"{jet_obs}{suffix}"
            if source not in columns:
                continue
            jet_obs_suff = "_".join(jet_obs.split("_")[1:])
            df = df.Define(
                f"vbfjet{vbfj_idx}_{jet_obs_suff}{suffix}",
                f"(HasVBF{suffix}) ? {source}.at(VBFJetIdx_{vbfj_idx}{suffix}): -1000.f;",
            )

    df = df.Define(
        f"m_jj{suffix}",
        f"if (HasVBF{suffix}) return static_cast<float>((vbfjet1_p4{suffix}+vbfjet2_p4{suffix}).M()); return -1000.f",
    )

    df = df.Define(
        f"delta_eta_jj{suffix}",
        f"if (HasVBF{suffix}) return static_cast<float>(abs(vbfjet1_p4{suffix}.Eta()-vbfjet2_p4{suffix}.Eta())); return -1000.f",
    )
    df = df.Define(
        f"vbfjet1_y{suffix}",
        f"if (HasVBF{suffix}) return static_cast<float>(vbfjet1_p4{suffix}.Rapidity()); return -1000.f; ",
    )
    df = df.Define(
        f"vbfjet2_y{suffix}",
        f"if (HasVBF{suffix}) return static_cast<float>(vbfjet2_p4{suffix}.Rapidity()); return -1000.f; ",
    )
    df = df.Define(
        f"delta_phi_jj{suffix}",
        f"if (HasVBF{suffix}) return static_cast<float>(ROOT::Math::VectorUtil::DeltaPhi(vbfjet1_p4{suffix},vbfjet2_p4{suffix})); return -1000.f;",
    )

    df = df.Define(
        f"pt_vbfj1j2{suffix}",
        f"if (HasVBF{suffix}) return static_cast<float>((vbfjet1_p4{suffix}+vbfjet2_p4{suffix}).Pt()); return -1000.f",
    )

    return df


def VBFJetMuonsObservablesDef(df):
    mu_suff = "ScaRe_FSR" if "mu1_p4_ScaRe_FSR" in df.GetColumnNames() else ""
    mu_suffix = f"_{mu_suff}" if mu_suff else ""
    mu1_p4 = f"mu1_p4{mu_suffix}"
    mu2_p4 = f"mu2_p4{mu_suffix}"
    pt_mumu_name = f"pt_mumu{mu_suffix}"
    eta_mumu_name = f"eta_mumu{mu_suffix}"
    y_mumu_name = f"y_mumu{mu_suffix}"

    for mu_idx in [1,2]:
        if f"mu{mu_idx}_p4" not in df.GetColumnNames():
            df = df.Define(f"mu{mu_idx}_p4", f"ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(mu{mu_idx}_pt, mu{mu_idx}_eta, mu{mu_idx}_phi, mu{mu_idx}_mass)")
    if "y_mumu" not in df.GetColumnNames():
        df = df.Define("y_mumu", "return (mu1_p4.Rapidity()+mu2_p4.Rapidity())/2.;")
    df = df.Define(
        "Zeppenfeld_Var",
        f"if (HasVBF) return static_cast<float>(({y_mumu_name} - 0.5*(vbfjet1_y+vbfjet2_y))/std::abs(vbfjet1_y - vbfjet2_y)); return -10000.f;",
    )
    df = df.Define(
        "pT_all_sum",
        f"if(HasVBF) return static_cast<float>(({mu1_p4}+{mu2_p4}+vbfjet1_p4+vbfjet2_p4).Pt()); return -10000.f;",
    )
    df = df.Define(
        "pT_single_sum",
        f"if(HasVBF) return static_cast<float>({mu1_p4}.Pt()+{mu2_p4}.Pt()+vbfjet1_p4.Pt()+vbfjet2_p4.Pt()); return -10000.f;",
    )
    df = df.Define(
        "R_pt",
        f"if(HasVBF) return static_cast<float>(pT_all_sum/pT_single_sum); return -10000.f;",
    )
    df = df.Define(
        "pT_jj_sum",
        f"if(HasVBF) return static_cast<float>((vbfjet1_p4+vbfjet2_p4).Pt()); return -10000.f;",
    )
    df = df.Define(
        "pT_jj_diff",
        f"if(HasVBF) return static_cast<float>((vbfjet1_p4-vbfjet2_p4).Pt()); return -10000.f;",
    )
    df = df.Define(
        "pT_mumu_sum",
        f"return static_cast<float>(({mu1_p4}+{mu2_p4}).Pt());",
    )
    df = df.Define(
        "pt_centrality",
        f"if(HasVBF) return static_cast<float>(( (pT_mumu_sum-0.5*(pT_jj_sum)) / pT_jj_diff)); return -10000.f;",
    )

    df = df.Define(
        "minDeltaPhi",
        f"if(HasVBF) return static_cast<float>(std::min(ROOT::Math::VectorUtil::DeltaPhi( ({mu1_p4}+{mu2_p4}), vbfjet1_p4), ROOT::Math::VectorUtil::DeltaPhi(({mu1_p4}+{mu2_p4}), vbfjet2_p4) ) )  ; return -10000.f;",
    )
    df = df.Define(
        "minDeltaEta",
        f"if(HasVBF) return static_cast<float>(std::min(std::abs({eta_mumu_name} - vbfjet1_eta),std::abs({eta_mumu_name} - vbfjet2_eta))) ; return -10000.f;",
    )
    df = df.Define(
        "minDeltaEtaSigned",
        f"if(HasVBF) return static_cast<float>(std::min(({eta_mumu_name} - vbfjet1_eta),({eta_mumu_name} - vbfjet2_eta))) ; return -10000.f;",
    )

    return df


def SoftJetCollectionCleaningInVBF(df):
    if "SoftActivityJet_pt" not in df.GetColumnNames():
        # Skims of the FlashSim Sep2026 ntuples lost the SoftActivityJet arrays
        # (ULong_t counters, see analysis/skim.py): only the cleaned soft-jet
        # observables are missing, the SoftActivityJetHT/Njets scalars are there.
        print("[WARNING] no SoftActivityJet_pt in the input: cleaned soft-jet observables not defined")
        return df
    if "SoftActivityJet_idx" not in df.GetColumnNames():
        # print("SoftActivityJet_idx not in df.GetColumnNames")
        df = df.Define(
            f"SoftActivityJet_idx", f"CreateIndexes(SoftActivityJet_pt.size())"
        )
    if f"SoftActivityJet_mass" not in df.GetColumnNames():
        df = df.Define(
            f"SoftActivityJet_mass",
            "RVecF SoftActivityJet_mass(SoftActivityJet_idx.size(),0.); return SoftActivityJet_mass;",
        )
    df = df.Define(
        f"SoftActivityJet_p4",
        f"GetP4(SoftActivityJet_pt, SoftActivityJet_eta, SoftActivityJet_phi, SoftActivityJet_mass, SoftActivityJet_idx)",
    )
    df = df.Define(
        f"SoftActivityJet_pt_gt0",
        f"SoftActivityJet_pt>0",
    )

    df = df.Define(
        "SoftJetActivity_NoOverlapWithMuonsAndVBFJets",
        f"""
        ROOT::VecOps::RVec<bool> out;
        out.reserve(SoftActivityJet_p4.size());

        for (size_t i = 0; i < SoftActivityJet_p4.size(); ++i)
        {{
            bool pass = ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], mu1_p4) > 0.4 &&
                    ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], mu2_p4) > 0.4 &&
                    ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], vbfjet1_p4) > 0.4 &&
                    ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], vbfjet2_p4) > 0.4;

            out.push_back(pass);
        }}

        return out;
        """,
    )

    df = df.Define(
        f"SoftJetCleanedActivity_pt",
        "v_ops::pt(SoftActivityJet_p4[SoftJetActivity_NoOverlapWithMuonsAndVBFJets])",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_eta",
        "v_ops::eta(SoftActivityJet_p4[SoftJetActivity_NoOverlapWithMuonsAndVBFJets])",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_N",
        "SoftActivityJet_p4[SoftJetActivity_NoOverlapWithMuonsAndVBFJets].size()",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_ptSum",
        "float sum=0.; for(size_t sj_idx=0; sj_idx<SoftJetCleanedActivity_pt.size();sj_idx++){{sum+=SoftJetCleanedActivity_pt[sj_idx];}} return sum;",
    )

    df = df.Define(
        "SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning",
        f"""
        ROOT::VecOps::RVec<bool> out;
        out.reserve(SoftJetCleanedActivity_eta.size());

        for (size_t i = 0; i < SoftJetCleanedActivity_eta.size(); ++i)
        {{
            bool pass = SoftJetCleanedActivity_eta[i] < std::max(vbfjet1_eta, vbfjet2_eta) && SoftJetCleanedActivity_eta[i] > std::min(vbfjet1_eta, vbfjet2_eta);

            out.push_back(pass);
        }}

        return out;
        """,
    )

    df = df.Define(
        f"SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_pt",
        "SoftJetCleanedActivity_pt[SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning]",
    )
    df = df.Define(
        f"SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_eta",
        "SoftJetCleanedActivity_eta[SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning]",
    )
    df = df.Define(
        f"SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_N",
        "SoftJetCleanedActivity_eta[SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning].size()",
    )
    df = df.Define(
        f"SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_ptSum",
        "float sum=0.; for(size_t sj_idx=0; sj_idx<SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_pt.size();sj_idx++){{sum+=SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_pt[sj_idx];}} return sum;",
    )
    return df


def _column_names(df):
    return {str(col) for col in df.GetColumnNames()}


def DefineDimuonMassResolution(df):
    """Define per-event dimuon mass resolution and ScaRe uncertainty terms.

    ``muN_pt_err`` is selected upstream from the beam-spot-constrained pT
    uncertainty when the BSC fit has chi2 < 30, and from the NanoAOD pT
    uncertainty otherwise.  Scale and smearing up/down shifts are systematic
    uncertainties, so they are propagated directly through the shifted
    dimuon masses instead of being folded into the detector resolution.
    """
    columns = _column_names(df)

    def define(name, expression):
        nonlocal df
        if name not in columns:
            df = df.Define(name, expression)
            columns.add(name)

    error_inputs = {"mu1_pt_err", "mu2_pt_err"}
    if error_inputs.issubset(columns):
        denominator_suffix = (
            "_noCorr"
            if {"mu1_pt_noCorr", "mu2_pt_noCorr"}.issubset(columns)
            else ""
        )
        denominator_columns = {
            f"mu1_pt{denominator_suffix}",
            f"mu2_pt{denominator_suffix}",
        }
        if denominator_columns.issubset(columns):
            for muon in (1, 2):
                pt = f"mu{muon}_pt{denominator_suffix}"
                define(
                    f"mu{muon}_pt_resolution_rel",
                    (
                        f"({pt} > 0.f && mu{muon}_pt_err >= 0.f) ? "
                        f"static_cast<float>(mu{muon}_pt_err / {pt}) : -1.f"
                    ),
                )
            define(
                "m_mumu_resolution",
                (
                    "(mu1_pt_resolution_rel >= 0.f && "
                    "mu2_pt_resolution_rel >= 0.f) ? "
                    "static_cast<float>(std::sqrt(0.5 * "
                    "(mu1_pt_resolution_rel * mu1_pt_resolution_rel + "
                    "mu2_pt_resolution_rel * mu2_pt_resolution_rel))) : -1.f"
                ),
            )
            # Explicit alias: m_mumu_resolution is dimensionless.
            define("m_mumu_resolution_detector", "m_mumu_resolution")
            if "m_mumu" in columns:
                define(
                    "m_mumu_resolution_abs",
                    (
                        "(m_mumu_resolution >= 0.f && m_mumu >= 0.f) ? "
                        "static_cast<float>(m_mumu_resolution * m_mumu) : -1.f"
                    ),
                )

    variation_pairs = {
        "scale": ("m_mumu_FSR_scale_up", "m_mumu_FSR_scale_down"),
        "res": ("m_mumu_FSR_res_up", "m_mumu_FSR_res_down"),
    }
    if "m_mumu" in columns:
        for component, (up, down) in variation_pairs.items():
            if {up, down}.issubset(columns):
                define(
                    f"m_mumu_resolution_{component}",
                    (
                        f"(m_mumu > 0.f) ? static_cast<float>"
                        f"(0.5 * std::abs({up} - {down}) / m_mumu) : -1.f"
                    ),
                )

    total_components = (
        "m_mumu_resolution_detector",
        "m_mumu_resolution_scale",
        "m_mumu_resolution_res",
    )
    if set(total_components).issubset(columns):
        valid = " && ".join(f"{name} >= 0.f" for name in total_components)
        quadrature = " + ".join(f"{name} * {name}" for name in total_components)
        define(
            "m_mumu_resolution_total",
            (
                f"({valid}) ? "
                f"static_cast<float>(std::sqrt({quadrature})) : -1.f"
            ),
        )
        if "m_mumu" in columns:
            define(
                "m_mumu_resolution_total_abs",
                (
                    "(m_mumu_resolution_total >= 0.f && m_mumu >= 0.f) ? "
                    "static_cast<float>(m_mumu_resolution_total * m_mumu) "
                    ": -1.f"
                ),
            )

    return df


def GetAllMuonsObservablesNew(df):
    df = df.Define("Ebeam", "13600.0/2")

    dimu_obs = {
        "pt_mumu": "{dimu}.Pt()",
        "m_mumu": "{dimu}.M()",
        "y_mumu": "{dimu}.Rapidity()",
        "eta_mumu": "{dimu}.Eta()",
        "phi_mumu": "{dimu}.Phi()",
        "dR_mumu": "ROOT::Math::VectorUtil::DeltaR({mu1p4}, {mu2p4})",
        "cosTheta_Phi_CS": "ComputeCosThetaPhiCS({mu1p4}, {mu2p4}, Ebeam)",
        "cosTheta_CS": "static_cast<float>(std::get<0>(cosTheta_Phi_CS{suff}))",
        "phi_CS": "static_cast<float>(std::get<1>(cosTheta_Phi_CS{suff}))",
    }
    for pt_suffix in ["", "_FSR_scale"
        # "_nano",
        # "_bsConstrainedPt",
        # "",  # should be same than bsc_scare
        # "_bsc_scare",
        # "_nano_scare",
        # "_FSR_nano_scare",
        # "_FSR_bsc_scare",
    ]:
        for mu_idx in [1, 2]:
            mu_pt_name = (
                f"mu{mu_idx}_pt{pt_suffix}"
                if pt_suffix != "_bsConstrainedPt"
                else f"mu{mu_idx}{pt_suffix}"
            )
            if (
                mu_pt_name in df.GetColumnNames()
                and f"mu{mu_idx}_p4{pt_suffix}" not in df.GetColumnNames()
            ):
                df = df.Define(
                    f"mu{mu_idx}_p4{pt_suffix}",
                    f"ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>({mu_pt_name},mu{mu_idx}_eta,mu{mu_idx}_phi,mu{mu_idx}_mass)",
                )
        p4_dimu = f"(mu1_p4{pt_suffix}+mu2_p4{pt_suffix})"
        p4_dimu_list = [f"mu1_p4{pt_suffix}", f"mu2_p4{pt_suffix}"]
        for obs, expr in dimu_obs.items():
            # if pt_suffix == "":
            #     continue
            if f"{obs}{pt_suffix}" in df.GetColumnNames():
                continue
            df = df.Define(
                f"{obs}{pt_suffix}",
                expr.format(
                    dimu=p4_dimu,
                    mu1p4=p4_dimu_list[0],
                    mu2p4=p4_dimu_list[1],
                    suff=pt_suffix,
                ),
            )

    columns = _column_names(df)
    # skim_v4 has no mu<i>_pt_noCorr: it keeps the two variants separate,
    # mu<i>_pt_FSR_noCorr (with FSR recovery) and mu<i>_pt_raw_noCorr (without).
    # The nominal m_mumu is the one with FSR -- its only systematic variants,
    # m_mumu_FSR_*, say so -- so the meaningful "without correction"
    # comparison is with the FSR variant: it isolates the effect of the
    # scale/resolution correction and not that of the FSR recovery. Before, the
    # block was simply skipped and m_mumu_noCorr stayed an empty histogram.
    no_corr_suffix = next(
        (
            suffix
            for suffix in ("_noCorr", "_FSR_noCorr", "_raw_noCorr")
            if all(f"mu{mu_idx}_pt{suffix}" in columns for mu_idx in [1, 2])
        ),
        None,
    )

    if no_corr_suffix is not None:
        for mu_idx in [1, 2]:
            p4_name = f"mu{mu_idx}_p4_noCorr"
            # eta and phi follow the same variant as the pT when they exist,
            # otherwise they stay the nominal ones: only the pT changes.
            angles = {
                axis: (
                    f"mu{mu_idx}_{axis}{no_corr_suffix}"
                    if f"mu{mu_idx}_{axis}{no_corr_suffix}" in columns
                    else f"mu{mu_idx}_{axis}"
                )
                for axis in ("eta", "phi")
            }
            if p4_name not in columns:
                df = df.Define(
                    p4_name,
                    (
                        f"ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>"
                        f"(mu{mu_idx}_pt{no_corr_suffix}, {angles['eta']}, "
                        f"{angles['phi']}, mu{mu_idx}_mass)"
                    ),
                )
                columns.add(p4_name)

        if "m_mumu_noCorr" not in columns:
            df = df.Define(
                "m_mumu_noCorr",
                "(mu1_p4_noCorr + mu2_p4_noCorr).M()",
            )
            columns.add("m_mumu_noCorr")

    df = DefineDimuonMassResolution(df)

    # for mu_idx in [1, 2]:
    #     df = df.Define(
    #         f"mu{mu_idx}_p4_noCorr",
    #         f"mu{mu_idx}_bsConstrainedChi2 < 30 ? mu{mu_idx}_p4_bsConstrainedPt : mu{mu_idx}_p4_nano",
    #     )
    #     df = df.Define(
    #         f"mu{mu_idx}_p4_ScaRe",
    #         f"mu{mu_idx}_bsConstrainedChi2 < 30 ? mu{mu_idx}_p4_bsc_scare : mu{mu_idx}_p4_nano_scare",
    #     )
    #     df = df.Define(
    #         f"mu{mu_idx}_p4_ScaRe_FSR",
    #         f"mu{mu_idx}_bsConstrainedChi2 < 30 ? mu{mu_idx}_p4_FSR_bsc_scare : mu{mu_idx}_p4_FSR_nano_scare",
    #     )
    # for newsuff in ["noCorr", "ScaRe", "ScaRe_FSR"]:
    #     df = df.Define(f"mu1_pt_{newsuff}", f"mu1_p4_{newsuff}.pt()")
    #     df = df.Define(f"mu2_pt_{newsuff}", f"mu2_p4_{newsuff}.pt()")
    #     p4_dimu_system = f"(mu1_p4_{newsuff}+mu2_p4_{newsuff})"
    #     p4_dimu_system_list = [f"mu1_p4_{newsuff}", f"mu2_p4_{newsuff}"]
    #     for obs, expr in dimu_obs.items():
    #         df = df.Define(
    #             f"{obs}_{newsuff}",
    #             expr.format(
    #                 dimu=p4_dimu_system,
    #                 mu1p4=p4_dimu_system_list[0],
    #                 mu2p4=p4_dimu_system_list[1],
    #                 suff=f"_{newsuff}",
    #             ),
    #         )
    #     df = df.Define(
    #         f"mu1_pt_rel_{newsuff}", f"mu1_p4_{newsuff}.pt()/m_mumu_{newsuff}"
    #     )
    #     df = df.Define(
    #         f"mu2_pt_rel_{newsuff}", f"mu2_p4_{newsuff}.pt()/m_mumu_{newsuff}"
    #     )

    # pt_variants = {
    #     "_nano": {
    #         "pt_err_template": "mu{0}_ptErr/mu{0}_pt",
    #         "pt_name_template": "mu{0}_pt_nano",
    #         "has_scare": False,
    #     },
    #     "_nano_scare": {
    #         "pt_err_template": "mu{0}_ptErr/mu{0}_pt",
    #         "pt_name_template": "mu{0}_pt_nano_scare",
    #         "has_scare": True,
    #         "base_p4_suffix": "_nano",
    #     },
    #     "_nano_scare_FSR": {
    #         "pt_err_template": "mu{0}_ptErr/mu{0}_pt",
    #         "pt_name_template": "mu{0}_pt_nano_scare_FSR",
    #         "has_scare": True,
    #         "base_p4_suffix": "_FSR_nano",
    #     },
    #     "_bsConstrainedPt": {
    #         "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
    #         "pt_name_template": "mu{0}_bsConstrainedPt",
    #         "has_scare": False,
    #     },
    #     "_bsc_scare": {
    #         "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
    #         "pt_name_template": "mu{0}_pt_bsc_scare",
    #         "has_scare": True,
    #         "base_p4_suffix": "_bsConstrainedPt",
    #     },
    #     "": {
    #         "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
    #         "pt_name_template": "mu{0}_pt_bsc_scare",
    #         "has_scare": True,
    #         "base_p4_suffix": "_bsConstrainedPt",
    #     },
    #     "_bsc_scare_FSR": {
    #         "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
    #         "pt_name_template": "mu{0}_pt_bsc_scare_FSR",
    #         "has_scare": True,
    #         "base_p4_suffix": "_FSR_bsConstrainedPt",
    #     },
    # }

    # for pt_suffix, pt_info in pt_variants.items():
    #     # Check if both muons have the required pT columns
    #     mu1_pt_name = pt_info["pt_name_template"].format(1)
    #     mu2_pt_name = pt_info["pt_name_template"].format(2)

    #     if (
    #         mu1_pt_name not in df.GetColumnNames()
    #         or mu2_pt_name not in df.GetColumnNames()
    #     ):
    #         continue

    #     # Calculate relative pT errors for each muon
    #     for mu_idx in [1, 2]:
    #         sigma_expr = pt_info["pt_err_template"].format(mu_idx)
    #         df = df.Define(f"sigma_mu{mu_idx}_pt_rel{pt_suffix}", sigma_expr)

    #     # Calculate m_mumu_resolution including ScaRe uncertainties
    #     # Base resolution from pT errors: Δm_μμ^rel = sqrt(1/2 * ((Δpt(u1)/pt(u1))^2 + (Δpt(u2)/pt(u2))^2))
    #     resolution_expr = f"sqrt(0.5*(pow(sigma_mu1_pt_rel{pt_suffix},2) + pow(sigma_mu2_pt_rel{pt_suffix},2)))"

    #     # Handle ScaRe uncertainties if present
    #     if pt_info.get("has_scare", False):
    #         base_p4_suffix = pt_info.get("base_p4_suffix", "")
    #         # Check if ScaRe delta columns exist from friend trees
    #         deltas_up_exist = all(
    #             f"mu{mu_idx}_p4{base_p4_suffix}_ScaReUp_delta" in df.GetColumnNames()
    #             for mu_idx in [1, 2]
    #         )
    #         deltas_down_exist = all(
    #             f"mu{mu_idx}_p4{base_p4_suffix}_ScaReDown_delta" in df.GetColumnNames()
    #             for mu_idx in [1, 2]
    #         )

    #         if deltas_up_exist and deltas_down_exist:
    #             # Use deltas from friend trees
    #             for mu_idx in [1, 2]:
    #                 mu_pt_name = pt_info["pt_name_template"].format(mu_idx)
    #                 scare_unc_name = f"sigma_mu{mu_idx}_pt_scare_rel{pt_suffix}"
    #                 delta_up = f"mu{mu_idx}_p4{base_p4_suffix}_ScaReUp_delta"
    #                 delta_down = f"mu{mu_idx}_p4{base_p4_suffix}_ScaReDown_delta"
    #                 # Average of absolute deltas, relative to nominal pT
    #                 df = df.Define(
    #                     scare_unc_name,
    #                     f"(abs({delta_up}) + abs({delta_down})) / 2.0 / {mu_pt_name}",
    #                 )

    #             # Add ScaRe uncertainties in quadrature to base resolution
    #             resolution_expr = f"sqrt(0.5*(pow(sigma_mu1_pt_rel{pt_suffix},2) + pow(sigma_mu2_pt_rel{pt_suffix},2)) + 0.5*(pow(sigma_mu1_pt_scare_rel{pt_suffix},2) + pow(sigma_mu2_pt_scare_rel{pt_suffix},2)))"

    #     resolution_name = (
    #         f"m_mumu_resolution{pt_suffix}" if pt_suffix != "" else "m_mumu_resolution"
    #     )
    #     df = df.Define(resolution_name, resolution_expr)

    return df


def ShiftedVBFMuonObservablesDef(df, jet_suffix="", muon_suffix="", out_suffix=None):
    """Derived VBF and dimuon observables for one jet or muon variation.

    VBFJetObservablesDef(suffix=...) already shifts the VBF-jet observables;
    this adds what depends on both jets and muons (Zeppenfeld_Var, R_pt,
    pt_centrality, minDeltaPhi, minDeltaEta[Signed]), the soft activity
    cleaned against the shifted VBF jets and, for a muon variation, the
    dimuon observables (pt_mumu, y_mumu, cosTheta_CS, ...).  Every column is
    named <observable><out_suffix>, with out_suffix = jet_suffix or
    muon_suffix by default, so the histogram maker and the DNN pick it up.
    With both suffixes empty and an explicit out_suffix it rebuilds the
    nominal columns, which is how the definitions are tested.
    """
    out = out_suffix if out_suffix is not None else (jet_suffix or muon_suffix)
    columns = _column_names(df)

    def define(name, expression):
        nonlocal df
        if name not in columns:
            df = df.Define(name, expression)
            columns.add(name)

    lorentz = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>"
    nominal_mu = "_ScaRe_FSR" if "mu1_p4_ScaRe_FSR" in columns else ""
    if muon_suffix:
        for mu_idx in (1, 2):
            angle = {axis: f"mu{mu_idx}_{axis}{muon_suffix}"
                     if f"mu{mu_idx}_{axis}{muon_suffix}" in columns else f"mu{mu_idx}_{axis}"
                     for axis in ("eta", "phi")}
            define(f"mu{mu_idx}_p4{out}",
                   f"{lorentz}(mu{mu_idx}_pt{muon_suffix}, {angle['eta']}, {angle['phi']}, mu{mu_idx}_mass)")
        mu1, mu2 = f"mu1_p4{out}", f"mu2_p4{out}"
        dimu = f"({mu1}+{mu2})"
        define("Ebeam", "13600.0/2")
        for name, expression in {
            "pt_mumu": f"{dimu}.Pt()",
            "y_mumu": f"{dimu}.Rapidity()",
            "eta_mumu": f"{dimu}.Eta()",
            "phi_mumu": f"{dimu}.Phi()",
            "dR_mumu": f"ROOT::Math::VectorUtil::DeltaR({mu1}, {mu2})",
            "cosTheta_Phi_CS": f"ComputeCosThetaPhiCS({mu1}, {mu2}, Ebeam)",
            "cosTheta_CS": f"static_cast<float>(std::get<0>(cosTheta_Phi_CS{out}))",
            "phi_CS": f"static_cast<float>(std::get<1>(cosTheta_Phi_CS{out}))",
        }.items():
            define(f"{name}{out}", expression)
        y_mumu, eta_mumu = f"y_mumu{out}", f"eta_mumu{out}"
    else:
        mu1, mu2 = f"mu1_p4{nominal_mu}", f"mu2_p4{nominal_mu}"
        y_mumu, eta_mumu = f"y_mumu{nominal_mu}", f"eta_mumu{nominal_mu}"

    j = jet_suffix
    has_vbf = f"HasVBF{j}"
    j1, j2 = f"vbfjet1_p4{j}", f"vbfjet2_p4{j}"
    define(f"Zeppenfeld_Var{out}",
           f"if ({has_vbf}) return static_cast<float>(({y_mumu} - 0.5*(vbfjet1_y{j}+vbfjet2_y{j}))"
           f"/std::abs(vbfjet1_y{j} - vbfjet2_y{j})); return -10000.f;")
    define(f"R_pt{out}",
           f"if ({has_vbf}) return static_cast<float>(({mu1}+{mu2}+{j1}+{j2}).Pt()"
           f"/({mu1}.Pt()+{mu2}.Pt()+{j1}.Pt()+{j2}.Pt())); return -10000.f;")
    define(f"pt_centrality{out}",
           f"if ({has_vbf}) return static_cast<float>((({mu1}+{mu2}).Pt()-0.5*({j1}+{j2}).Pt())"
           f"/({j1}-{j2}).Pt()); return -10000.f;")
    define(f"minDeltaPhi{out}",
           f"if ({has_vbf}) return static_cast<float>(std::min("
           f"ROOT::Math::VectorUtil::DeltaPhi(({mu1}+{mu2}), {j1}), "
           f"ROOT::Math::VectorUtil::DeltaPhi(({mu1}+{mu2}), {j2}))); return -10000.f;")
    define(f"minDeltaEta{out}",
           f"if ({has_vbf}) return static_cast<float>(std::min(std::abs({eta_mumu} - vbfjet1_eta{j}),"
           f"std::abs({eta_mumu} - vbfjet2_eta{j}))); return -10000.f;")
    define(f"minDeltaEtaSigned{out}",
           f"if ({has_vbf}) return static_cast<float>(std::min(({eta_mumu} - vbfjet1_eta{j}),"
           f"({eta_mumu} - vbfjet2_eta{j}))); return -10000.f;")

    # Soft activity: the cleaning against the VBF jets moves with them.  A
    # muon variation changes only the muon pT, not the cleaning cones.
    if jet_suffix or out_suffix is not None:
        if "SoftActivityJet_p4" in columns:
            define(f"SoftJetActivity_NoOverlapWithMuonsAndVBFJets{out}", f"""
                ROOT::VecOps::RVec<bool> keep(SoftActivityJet_p4.size());
                for (size_t i = 0; i < SoftActivityJet_p4.size(); ++i)
                    keep[i] = ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], mu1_p4) > 0.4 &&
                              ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], mu2_p4) > 0.4 &&
                              ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], {j1}) > 0.4 &&
                              ROOT::Math::VectorUtil::DeltaR(SoftActivityJet_p4[i], {j2}) > 0.4;
                return keep;""")
            mask = f"SoftJetActivity_NoOverlapWithMuonsAndVBFJets{out}"
            define(f"SoftJetCleanedActivity_N{out}", f"SoftActivityJet_p4[{mask}].size()")
            define(f"SoftJetCleanedActivity_ptSum{out}",
                   f"float sum = 0.; for (const auto& p4 : SoftActivityJet_p4[{mask}]) sum += p4.Pt(); return sum;")
    return df
