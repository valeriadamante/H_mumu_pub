# TAMU sync tuple (`update_after_Hmm_15Oct`): variables and how they are filled

Written 2026-10-19 from `Sync_tuple_MC_Types.txt` (635 columns of `TAMU-GluGluHto2Mu.root`, tree `Events`) and from the code that fills them:

- `sync/python/special_sync_skim_update.py` (NanoAOD -> tuple; muons, dimuon, electrons, jets, VBF pairs, weights, blinding)
- `config/Run3_2024/selections_sync_ggfvbf_update_20261016.yaml` (all regions, flags and categories)
- `sync/python/prepare_tuple_forSync.py` (`prepare`, `selected_columns`: leading/subleading/third/VBF jet columns, dijet columns)
- `analysis/muons.py` (muon pT/p4 variants), `corrections/general.py` (`define_base_weights`)

Items marked **(not verified)** are inferred from names or from code that is not in this repository's sync scripts.
Missing values: the per-jet columns (`leadingjet_*`, `vbfjet1_*`, ...) are `-1000` when the jet does not exist; `vbf_*_jet*_*`, `mjj_*`, `deta_jj_*` and `dijet_*` use `-999`/`-1000` as noted below.
The columns are the same for data and MC except for the MC-only `weight_*`, `genWeight` and `mu*_GenMatched`/`mu*_genPart*` columns.

> **Note on `mu{i}_pt`.** In this tuple `mu{i}_pt` and `mu{i}_eta` are the plain **NanoAOD** values (`Muon_pt[idx]`, `Muon_eta[idx]`): they are set back to NanoAOD after being temporarily redefined to the FSR-corrected values for the MC scale factors (`special_sync_skim_update.py:459-477`).
> The corrected muon kinematics (BSC/NanoAOD switch + ScaRe + FSR, the "production" pT used to order mu1/mu2 and to build `m_mumu`) are `mu{i}_pt_FSR_corr`, `mu{i}_eta_FSR_corr`, `mu{i}_phi_FSR_corr`.
> Always compare corrected muon pT with `mu{i}_pt_FSR_corr` (as `compare_run3_sync_groups.py:117-118` does), and the NanoAOD pT with `mu{i}_pt` (= `mu{i}_pt_nano`). The name is misleading; renaming is pending.

## 0. Selection of the tuple (what each row has already passed)

Events are kept if `skim_selection` is true: `base_sel_muons && mass_50_200`, where

- `base_sel_muons = good_event && muons_OS_presel_trg && mu_pt_sel && bothMuons_ML && exactly_2_good_muons && electron_veto`
- `good_event = pass_golden_json && pass_met_filters` (for MC `pass_golden_json` is always true)
- `mass_50_200 = 50 < m_mumu_raw < 200`
- Data only: events with any of `m_mumu_raw`, `m_mumu_nano`, `m_mumu_nano_corr`, `m_mumu_FSR_corr`, `m_mumu` in 115-135 GeV are removed (blinding). MC is not blinded and is not lumi-split.
- The b-veto is **not** a filter: it is the flag `SelectedJetTagSel`/`jets_btag_clean_sel`, applied in `baseline` and therefore in `VBF`, `ggF`, ...

## 1. Event identification and bookkeeping

| Column | Type | How it is filled |
|---|---|---|
| `run`, `luminosityBlock`, `event` | UInt/ULong | copied from NanoAOD |
| `sample` | string | dataset name given to the job (e.g. `GluGluHto2Mu`) |
| `sample_id` | Int | `abs(crc32(dataset name))` |
| `is_data` | Bool | from `samples.yaml` (`is_data`) |
| `period` | Int | era enum of the config (`Period::Run3_2024`) |
| `MET_pt`, `MET_phi` | Float | `PuppiMET_pt/phi` |
| `MET_collection` | string | constant `"PuppiMET"` |
| `pass_golden_json` | Bool | data: golden-JSON lumi mask on (run, lumi); MC: `true` |
| `pass_met_filters` | Bool | AND of the `MET_flags` list of `maincfg.yaml`, plus the bad-MET-run filter if configured |
| `Flag_*` (about 28 columns) | Bool | the NanoAOD MET/noise filter flags, copied unchanged |
| `HLT_IsoMu24` | Bool | NanoAOD trigger bit |
| `Event_HasTriggerMatching_singleMu` | Bool | `HLT_IsoMu24 && any muon matched (dR < 0.4, NanoAOD pT and direction) to a trigger object passing the online cut`; the offline cut applied to the matched muon is `Muon_pt > 26` (canonical matching) |
| `sync_preselection` | Bool | `m_mumu_raw` in 50-200 |

## 2. Muons (`mu1_*`, `mu2_*`)

**Which muons.** Good muon = `Muon_pt_raw_noCorr > 20 && |Muon_eta| < 2.4 && mediumId && pfRelIso04_all < 0.25`. Exactly 2 good muons are required. The "exactly 2 loose muons" filter exists in the code but is not applied in the update (inclusive) variant; `n_looseMuons` is only stored. `mu1`/`mu2` are the two good muons ordered by **`Muon_pt_FSR_corr`** (BSC/NanoAOD pT + ScaRe + FSR), highest first. `mu{i}_idx` is the index into the NanoAOD `Muon` collection.

**pT variants.** The BSC choice is made per muon: `Muon_bsConstrainedChi2 < 30` → use the beam-spot-constrained (BSC) quantity, otherwise the NanoAOD one.

| Column | How it is filled |
|---|---|
| `mu{i}_pt`, `mu{i}_eta`, `mu{i}_phi`, `mu{i}_mass` | **NanoAOD** `Muon_pt/eta/phi/mass` of the muon (no ScaRe, no FSR, no BSC). In the MC weight step they are temporarily redefined to the FSR-corrected values to evaluate the scale factors, then reset to the NanoAOD ones. Note: the docstring of `prepare_tuple_forSync.py` (`CONVENTIONS`) still says `mu_pt` is the corrected one; this is not the case in the update tuple |
| `mu{i}_pt_nano` | same as `mu{i}_pt` (`Muon_pt`) |
| `mu{i}_pt_nano_corr` | `Muon_pt_nano_corr`: NanoAOD pT with ScaRe |
| `mu{i}_pt_nano_corr_FSR` | `Muon_pt_nano_corr_FSR`: NanoAOD pT with ScaRe and FSR recovery |
| `mu{i}_pt_raw_noCorr` | BSC if `chi2 < 30`, else NanoAOD; no corrections |
| `mu{i}_pt_raw_scale` | same choice, ScaRe scale correction only **(not verified)** |
| `mu{i}_pt_raw_corr` | same choice, ScaRe scale + resolution smearing (smearing in MC) |
| `mu{i}_pt_FSR_noCorr` | same choice, FSR photon added, no ScaRe |
| `mu{i}_pt_FSR_scale` | same choice, FSR + ScaRe scale only **(not verified)** |
| `mu{i}_pt_FSR_corr` | same choice, ScaRe + FSR: **the production pT** |
| `mu{i}_eta/phi_{raw_noCorr,raw_corr,raw_scale,FSR_noCorr,FSR_corr,FSR_scale}` | eta/phi of the p4 `Muon_p4_<stage>` of that muon (with FSR, the photon-added direction); Double |
| `mu{i}_pt_err` | `Muon_pt_err`: BSC pT error if `chi2 < 30`, else NanoAOD `Muon_ptErr` |
| `mu1_pt_g26`, `mu2_pt_g20` | flags `mu1_pt_FSR_corr > 26`, `mu2_pt_FSR_corr > 20`; stored for comparison only, not used by any category |
| `mu{i}_bsConstrainedChi2` | NanoAOD `Muon_bsConstrainedChi2` |
| `mu{i}_charge`, `mu{i}_pdgId` | NanoAOD |
| `mu{i}_dxy`, `_dxyErr`, `_dz`, `_dzErr`, `_dxybs`, `_dxybsErr`, `_IPx/_IPy/_IPz`, `_VXBS_Cov00/03/33` | NanoAOD muon impact-parameter and beam-spot vertex quantities, copied from `Muon_<name>[idx]` |
| `mu{i}_looseId`, `_mediumId`, `_tightId`, `_highPtId`, `_pfIsoId`, `_miniIsoId`, `_pfRelIso04_all`, `_miniPFRelIso_all`, `_tkRelIso`, `_nTrackerLayers` | NanoAOD ID/isolation variables |
| `mu{i}_newIso` | `(pfRelIso04_all * Muon_pt - pT of the associated FSR photon inside the isolation cone) / Muon_pt` (FSR photon inside dR < 0.4, see below). **Diagnostic: not used** in the update selection, which uses the raw `pfRelIso04_all` |
| `mu{i}_fsrPhotonIdx` | NanoAOD index of the FSR photon associated with the muon |
| `mu{i}_fsrPhoton_ptIso` | pT of that photon if it passes the FSR criteria and `0.0001 < dR(mu, photon) < 0.4`, else 0 |
| `mu{i}_fsrPhoton_ptRecovered` | same with the full recovery acceptance, `dR < 0.5` |
| `mu{i}_fsrPhoton_dR` | dR(muon, associated photon), -1 if none. The photon criteria are: no electron match, `pT_photon/pT_mu < 0.4`, `dROverEt2 < 0.012`, `relIso03/pT_mu < 1.8` |
| `mu{i}_svIdx` | NanoAOD secondary-vertex index |
| `mu{i}_HasTriggerMatching_singleMu` | muon matched (dR < 0.4) to the single-muon trigger object |
| `mu{i}_GenMatched` | MC: `Muon_genPartFlav` is 1 (prompt) or 15 (from tau); data: false |
| `mu{i}_genPartFlav`, `mu{i}_genPartIdx` | NanoAOD gen matching (MC) |
| `mu{i}_ML` | `mediumId && pfRelIso04_all < 0.25` |
| `bothMuons_ML` | `mu1_ML && mu2_ML` |
| `n_goodMuons`, `n_looseMuons` | number of good muons (above) and of loose muons (`pT_raw_noCorr > 20, |eta| < 2.4, looseId, pfRelIso04_all < 0.25`) |
| `mu_pt_sel` | `mu1_pt_raw_noCorr > 20 && mu2_pt_raw_noCorr > 20` |
| `muon_pre_sel_check` | raw pT > 20 and `|eta_raw_noCorr| < 2.4` for both muons |
| `muons_OS` | `mu1_charge * mu2_charge < 0` |
| `muons_trigger` | `Event_HasTriggerMatching_singleMu && HLT_IsoMu24 && (mu1 matched || mu2 matched)` |
| `chi2_sel` | both `bsConstrainedChi2 < 30` (stored only, not used) |
| `exactly_2_good_muons` | `n_goodMuons == 2` |

## 3. Dimuon system

All masses are for the pair (mu1, mu2). `dimuon_*` are aliases of `m_mumu`, `pt_mumu`, `eta_mumu`, `phi_mumu`.

| Column | Four-vectors used |
|---|---|
| `m_mumu`, `dimuon_mass` | `Muon_p4_FSR_corr`: BSC/NanoAOD switch + ScaRe + FSR (**same as `m_mumu_FSR_corr`**) |
| `pt_mumu`, `eta_mumu`, `phi_mumu` (`dimuon_pt/eta/phi`) | same four-vector as `m_mumu` |
| `m_mumu_FSR_corr` | `(mu1_p4_FSR_corr + mu2_p4_FSR_corr).M()` |
| `m_mumu_raw` | `Muon_p4_raw_noCorr`: BSC/NanoAOD switch, no corrections. Used for the 50-200 window and for the blinding |
| `m_mumu_nano` | NanoAOD pT/eta/phi/mass, no corrections |
| `m_mumu_nano_corr` | NanoAOD pT with ScaRe only |
| `m_mumu_nano_corr_FSR` | NanoAOD pT with ScaRe and FSR |
| `pt_mumu_nano` | pT of the NanoAOD-only pair |

## 4. Mass regions (flags, on `dimuon_mass` = `m_mumu`)

| Flag | Definition |
|---|---|
| `Z_sideband` | 70 < m < 110 |
| `Signal_Fit` | 115 < m < 135 (empty in data after the blinding) |
| `H_sideband` | 110 < m < 115 or 135 < m < 150 |
| `Signal_ext` | 110 < m < 150 |
| `mass_inclusive` | 60 < m < 150 |
| `Z_CR_H_SB` | 76 < m < 115 or 130 < m < 150 |
| `mass_50_200` | 50 < `m_mumu_raw` < 200 (raw mass) |

## 5. Electrons

| Column | How it is filled |
|---|---|
| `n_goodElectrons` | `Electron_pt > 10 && |eta| < 2.5 && mvaIso_WP90` |
| `n_goodElectrons_pt20` | same with `pT > 20` |
| `n_goodElectrons_WP80` | `pT > 10, |eta| < 2.5, mvaIso_WP80` |
| `n_goodElectrons_dRmu04` | good electrons (WP90, pT > 10) farther than 0.4 from every good muon |
| `electron_veto` | `n_goodElectrons_pt20 == 0` (**the veto of the update**) |

## 6. Jets

Jet collections come from the framework (JEC reapplied, JER in MC). The definitions are those of the update: pT > 30, |eta| < 4.7, tight jet ID (`Jet_passJetIdTight`), outside the veto maps, outside the 2024 horn region, no overlap (dR > 0.4) with mu1/mu2 (`special_sync::CleanJets`).

| Column family | How it is filled |
|---|---|
| `Jet_*` (vectors, all jets) | the full jet collection with the framework columns (pT/mass corrected, `_nocorr`, `_pt_raw`, b-tag scores, ID, energy fractions, `hadronFlavour`, `partonFlavour`, `genJetIdx`, ...) |
| `Jet_pt`, `Jet_mass` | JEC reapplied (+ JER in MC) |
| `Jet_pt_nocorr`, `Jet_mass_nocorr` | NanoAOD jet pT/mass before our JEC reapplication and JER |
| `Jet_pt_raw` | `Jet_pt_nocorr * (1 - Jet_rawFactor)` (no JEC) |
| `Jet_btag_loose`, `Jet_btag_medium` | UParTAK4 b-tag working points, only |eta| < 2.5 |
| `Jet_vetoMap`, `Jet_vetoMapLooseRegion_presel`, `Jet_isInsideVetoRegion` | jet in the JME veto map |
| `Jet_IsInsideHorn`, `Jet_IsOutsideHorn` | horn region of the config (2024: 2.5 <= |eta| < 3.0 and pT <= 50 GeV) |
| `Jet_preSel` | `pT > 30 && |eta| < 4.7 && passJetIdTight` |
| `Jet_NoOverlapWithMuons` | `Jet_preSel && !vetoMap && outside horn && dR(jet, mu1/mu2) > 0.4` (electrons are not used for the cleaning) |
| `Jet_passJetIdTight`, `Jet_passJetIdTightLepVeto` | tight jet ID recomputed from the energy fractions |
| `Jet_idx`, `SelectedJet_idx` | index in the full jet collection |
| `SelectedJet_*` (vectors) | the same quantities for the jets with `Jet_NoOverlapWithMuons`, ordered by corrected pT |
| `N_SelectedJets` | number of selected jets |
| `nJet` | number of jets in the `Jet` collection (not the selected ones) **(not verified)** |
| `leadingjet_*`, `subleadingjet_*`, `thirdleadingjet_*` | element 0, 1, 2 of every `SelectedJet_*` vector, as Float; `-1000` if the event has fewer jets |
| `SelectedJetTagSel` (= `jets_btag_clean_sel`), `JetTagSel` | b-veto flag: no medium and fewer than two loose UParTAK4 b-tagged selected jets (|eta| < 2.5). `JetTagSel` is the version on all jets **(not verified)** |
| `has_good_jets` | `N_SelectedJets >= 0` (always true; placeholder) |
| `dijet_ls_mass/pt/eta/phi` | two leading selected jets; `-1000` if fewer than 2 jets |
| `mjj_j0j1`, `deta_jj_j0j1` | m_jj and |deta| of the two leading selected jets whatever the VBF cuts (`-999` if fewer than 2) |

## 7. VBF pairs and categories

Two VBF definitions are stored side by side.

| Quantity | Leading-jet definition (default) | Max-m_jj definition |
|---|---|---|
| Pair | the two leading selected jets (index 0, 1) | among all pairs (i < j) the one with the highest m_jj |
| Cuts | pT1 > 35, pT2 > 30, `|deta| > 2.5`, `m_jj > 400` | same cuts |
| Flag | `HasVBF` = `HasVBF_leading` = `VBF_def` | `HasVBF_maxmjj` = `VBF_def_maxmjj` |
| Jet indices | `VBFJetIdx_1/2` = `VBFJetIdx_leading_1/2` (-1000 if none) | `VBFJetIdx_maxmjj_1/2` |
| Kinematics | `vbf_leading_jet{1,2}_{pt,eta,phi,mass}`, `mjj_leading`, `deta_jj_leading` (-999 if none) | `vbf_maxmjj_jet{1,2}_*`, `mjj_maxmjj`, `deta_jj_maxmjj` |
| Categories | `VBF`, `ggF`, `ggF_0J/1J/2J/ge2J`, `VBF_ge2J` | `VBF_maxmjj`, `ggF_maxmjj`, `VBF_maxmjj_ge2J` |

- `VBF_maxmjj_same_pair_as_leading`: both flags true and the max-m_jj pair is (0, 1).
- `vbfjet1_*`, `vbfjet2_*`: every `SelectedJet_*` quantity for the jets `VBFJetIdx_1/2` of the **leading** definition (`-1000` if no pair).
- `dijet_vbf_mass/pt/eta/phi`: four-vector sum of the two jets `VBFJetIdx_1/2`; `-1000` if no pair.
- The categories are exclusive between VBF and ggF within `baseline`:
  - `baseline = base_sel_muons && SelectedJetTagSel` (this is where the b-veto enters)
  - `VBF = baseline && VBF_def`
  - `ggF = baseline && !VBF_def`
  - `ggF_0J/1J/2J/ge2J = ggF && N_SelectedJets == 0 / 1 / 2 / >= 2`
  - `VBF_ge2J = VBF && N_SelectedJets >= 2`

## 8. MC weights (MC only)

| Column | How it is filled |
|---|---|
| `genWeight`, `weight_generator` | NanoAOD generator weight (the second is a copy) |
| `weight_gen_sign` | `copysign(1, genWeight)` |
| `weight_lumi` | integrated luminosity of the era from `maincfg.yaml` (the 2024 value, as a constant float) |
| `weight_xs` | cross section from `crossSections13p6TeV.yaml` (for DY-like processes the DY cross-section processor) |
| `weight_pu`, `weight_pileup` | pileup weight from `Pileup_nTrueInt` (correctionlib); the second is a copy |
| `weight_mu{i}_<ID>_<Iso>_Central` (about 18 per muon) | muon ID/iso/tracking scale factors from the muon POG `muon_Z.json.gz`, central value, evaluated with the FSR-corrected pT and eta of the muon |
| `weight_TrgSF_singleMu_IsoMu24Central`, `weight_trigger` | IsoMu24 trigger SF (the second is a copy) |
| `weight_muon_idiso` | product over the two muons of ID/iso SF and tracking SF: for pT < 30 GeV the Tight ID + TightPFIso (+ TightID_Trk) factors, otherwise MediumID + LoosePFIso (+ MediumID_Trk) |

**Normalisation (to be confirmed).** The tuple does not store the factor `1/sum(genWeight)` of the whole dataset (the job reports do not contain it either), so the normalised event weight is
`genWeight * weight_lumi * weight_xs * weight_pu * weight_muon_idiso * weight_trigger / sum(genWeight)` with `sum(genWeight)` taken over the **full** dataset (all jobs, before any selection), not over the tuple. The b-tag weight is not in the tuple. Check this combination against `compare_run3_sync_groups.py` before using it.

## 9. Differences to keep in mind when comparing with other groups

- The corrected muon pT of TAMU is `mu{i}_pt_FSR_corr` (and the corrected dimuon mass is `m_mumu` = `m_mumu_FSR_corr`); `mu{i}_pt` is the NanoAOD pT.
- Muons are ordered by `Muon_pt_FSR_corr`, so mu1 is the leading muon after ScaRe + FSR.
- `leadingjet_*` is defined with at least 1 selected jet, `subleadingjet_*` with at least 2.
- The data tuple has no events in 115-135 GeV for any of the mass definitions.
