# Tupla loose per il sync con Pisa (Run3_2024)

Aggiornato: 9 ottobre 2026.

Tupla con selezione larga, senza andare subito nelle regioni ggF/VBF: tutti gli
eventi dello skim_v4 con 50 < m_μμ < 200 GeV, con le definizioni di regioni,
categorie e tagli sui muoni salvate come flag. Così una differenza con Pisa si
può ricondurre al singolo taglio senza rifare lo skim.

Codice: `sync/python/prepare_tuple_forSync.py` (tupla, cutflow, metadata) e
`sync/python/add_raw_jet_pt.py` (pT dei jet non corretto, lanciato in automatico).
Definizioni dei flag: `config/Run3_2024/selections_sync.yaml`.

## Produzione

Da lanciare in `tmux`/`screen`: il data richiede un paio d'ore.

```bash
source env.sh
SKIM=$HMUMU_SKIM/Run3_2024
OUT=$HMUMU_SYNC/loose_tuples_20261009

# data Run2024H (Muon0 + Muon1 in un'unica tupla)
python3 sync/python/prepare_tuple_forSync.py $SKIM/Muon0_Run2024H $SKIM/Muon1_Run2024H \
  -o $OUT/sync_tuple_loose_Run2024H.root --threads 12

# MC: un dataset per tupla (la normalizzazione è quella del dataset intero)
for ds in VBFHto2Mu_M125_powheg GluGluHto2Mu; do
  python3 sync/python/prepare_tuple_forSync.py $SKIM/$ds --mc \
    -o $OUT/sync_tuple_loose_$ds.root --threads 12
done
```

Per altri campioni basta cambiare il nome del dataset: è il nome della cartella
dello skim, e deve esistere anche in `config/Run3_2024/samples_withfiles.yaml`,
da cui si prendono i NanoAOD per il pT raw dei jet.

Opzioni utili:

| opzione | effetto |
|---|---|
| `--selection EXPR` | filtro della tupla (default `mass_50_200`; `1` = nessun filtro) |
| `--no-raw-jet-pt` | salta la lettura dei NanoAOD (niente `*_pt_raw` dei jet) |
| `--nano-workers N` | letture NanoAOD in parallelo (default 24) |
| `--no-framework` | niente `prepare_rdf`: niente pesi né osservabili del framework, più veloce |
| `--overwrite` | sovrascrive una tupla esistente |
| `--cutflow-key KEY` | lista di step del cutflow in `selections_sync.yaml` |

### Tempi indicativi

| passo | VBFH (11.2 M eventi NanoAOD, 3.4 M nello skim) | data Run2024H (188 M NanoAOD, 4.3 M nello skim) |
|---|---|---|
| tupla da skim | ~50 s per un file skim di prova (soprattutto overhead); dataset intero non misurato | non misurato, 60 file skim |
| lettura NanoAOD (`Jet_rawFactor`) | ~5 min, 190 file, 24 worker | ~1.5 h stimate, 245 file |

La lettura NanoAOD è salvata in `<output>_nano_cache/`: se si interrompe, rilanciare
lo stesso comando con `--overwrite` riprende dai file già letti. Quando la tupla è
definitiva la cartella si può cancellare.

Il pT raw si può aggiungere anche dopo, a una tupla già fatta:

```bash
python3 sync/python/add_raw_jet_pt.py $OUT/sync_tuple_loose_Run2024H.root \
  --dataset Muon0_Run2024H --dataset Muon1_Run2024H
```

## File prodotti

Per ogni `-o NAME.root`:

| file | contenuto |
|---|---|
| `NAME.root` | albero `Events` + `TNamed metadata` (lo stesso JSON qui sotto) |
| `NAME_metadata.json` | dataset, lumi e sezioni d'urto, correzioni, selezioni, convenzioni |
| `NAME_cutflow.csv` | cutflow dai `report_*.json` dello skim (dal NanoAOD) fino a `mass_inclusive` |
| `NAME_categories_regions.csv` | eventi per ogni regione/categoria salvata |

## Contenuto della tupla

Gli indici dei muoni sono `mu1` (leading) e `mu2`. Il valore mancante è −1000.
Sono salvate solo le quantità nominali: le varianti JES/JER e FSR scale/res e i
pesi up/down sono esclusi.

**Evento**: `run`, `luminosityBlock`, `event`, `HLT_IsoMu24`,
`Event_HasTriggerMatching_singleMu`, `is_data`, `period`, `OrthogonalEraTag`.

**Flag di regione** (massa con muoni ScaRe+FSR): `mass_50_200`, `mass_inclusive`
(60–150), `Z_sideband` (70–110), `Z_CR_H_SB`, `H_sideband`, `Signal_Fit` (115–135),
`Signal_ext` (110–150).

**Flag di categoria**: `base_sel_muons`, `baseline`, `VBF_def`, `VBF`, `VBF_ge2J`,
`ggF`, `ggF_0J`, `ggF_1J`, `ggF_2J`, `ggF_ge2J`.

**Flag sui muoni e sui jet**: `muon_pre_sel_check`, `muons_OS`, `muons_trigger`,
`mu1_rawpt_g26`, `mu2_rawpt_g20`, `mu_rawpt_sel` (pT raw), `mu1_pt_g26`, `mu2_pt_g20`,
`mu2_pt_g26` (pT corretto, solo per confronto), `mu1_ML`, `mu2_ML`, `bothMuons_ML`,
`chi2_sel`, `has_good_jets`, `jets_btag_clean_sel`, `SelectedJetTagSel`, `JetTagSel`.

Tagli sui muoni del sync (come Pisa): leading > 26 GeV e subleading > 20 GeV sul
**pT raw** `mu*_pt_raw_noCorr` (`mu_rawpt_sel`). Sono dentro `base_sel_muons`, quindi
anche in `baseline`, `VBF` e `ggF*`. La `baseline` del sync non ha il requisito
tight-tight sotto 30 GeV. La versione con gli stessi tagli sul pT corretto è
`mu1_pt_g26 && mu2_pt_g20`. Fino all'8 ottobre la baseline del sync era 26/26 sul
pT corretto.

Due limiti dello skim_v4, che non si possono recuperare da questa tupla:

- il "raw" dello skim è il pT senza ScaRe né FSR, ma con la scelta BSC
  (`bsConstrainedPt` se χ² < 30). Il pT NanoAOD puro `Muon_pt`, quello che usa
  Pisa, non è salvato;
- lo skim richiede che il muone associato al trigger abbia pT raw > 26 GeV.
  Pisa chiede solo un muone leading > 26 e almeno un muone associato. Nel primo
  confronto Run2024H questo spiegava 103 eventi solo-Pisa
  (`LESSONS_LEARNED_PISA.md`, punto 5).

**Muoni** (`mu1_*`, `mu2_*`):

| branch | definizione |
|---|---|
| `mu*_pt_raw_noCorr` | pT NanoAOD senza ScaRe né FSR: `bsConstrainedPt` se `bsConstrainedChi2 < 30`, altrimenti `Muon_pt` |
| `mu*_pt_raw_scale` | stessa scelta + scala ScaRe |
| `mu*_pt_raw_corr` | stessa scelta + scala ScaRe + smearing (MC) |
| `mu*_pt` (= `mu*_pt_FSR_corr`) | nominale: ScaRe + recupero FSR; `mu*_eta/phi/mass` dallo stesso quadrivettore |
| `mu*_eta_raw_noCorr`, `mu*_phi_raw_noCorr` | angoli NanoAOD |
| `mu*_mediumId`, `tightId`, `looseId`, `highPtId`, `pfIsoId`, `pfRelIso04_all`, `miniIsoId`, `tkRelIso` | ID e isolamento |
| `mu*_charge`, `dxy`, `dz`, `bsConstrainedChi2`, `pt_err`, `fsrPhotonIdx` | altro |
| `mu*_HasTriggerMatching_singleMu` | muone associato all'oggetto di trigger |
| `mu*_GenMatched`, `genPartFlav` | solo MC |

**Dimuone**: `m_mumu`/`dimuon_mass`, `pt_mumu`/`dimuon_pt`, `eta_mumu`, `phi_mumu`,
`y_mumu`, `dR_mumu` (muoni ScaRe+FSR); `m_mumu_noCorr`, `pt_mumu_noCorr` (pT
`raw_noCorr` e angoli NanoAOD).

**Jet**, tre versioni del pT:

| branch | definizione |
|---|---|
| `Jet_pt`, `SelectedJet_pt` | nostra JEC riapplicata (+ JER nel MC) |
| `Jet_pt_nocorr`, `SelectedJet_pt_nocorr` | pT NanoAOD (JEC di Nano), prima della nostra JEC |
| `Jet_pt_raw`, `SelectedJet_pt_raw` | pT raw, senza JEC né JER: `Jet_pt_nocorr × (1 − Jet_rawFactor)` |

- `Jet_*`: tutti i jet nell'ordine NanoAOD. `SelectedJet_*`: i jet selezionati, ordinati per pT corretto; `SelectedJet_idx` è l'indice in `Jet_*`.
- Qualità: `passJetIdTight`, `passJetIdTightLepVeto` (ID ricalcolato dalle frazioni di energia, ricetta NanoAOD v15), `vetoMap`, `isInsideVetoRegion`, `IsInsideHorn`/`IsOutsideHorn`, `NoOverlapWithMuons`, `preSel`, `puIdDisc`.
- b-tag: `btagUParTAK4B` (tagger usato) e gli altri score; `btag_loose`, `btag_medium` (WP UParTAK4, solo |η| < 2.5).
- Jet appiattiti: `leadingjet_*`, `subleadingjet_*`, `thirdjet_*`, `thirdleadingjet_*` (pT, η, φ, massa, b-tag, …) e `*_pt_raw`; `N_SelectedJets`, `nJet`.
- Dijet: `dijet_ls_{mass,pt,eta,phi}` e `m_jj_ls`, `delta_eta_jj_ls` (due jet leading); `HasVBF`, `VBFJetIdx_1/2`, `vbfjet1_*`, `vbfjet2_*`, `dijet_vbf_{mass,pt,eta,phi}`, `m_jj`, `delta_eta_jj`, `delta_phi_jj`, `pt_vbfj1j2` (coppia VBF).
- Solo MC: `hadronFlavour`, `partonFlavour`, `genJetIdx`.

**Pesi (solo MC)**:

```
weight__Central = genWeight × weight_lumi × weight_xs × weight_pu
                × SF ID × SF iso (mu1, mu2) × weight_TrgSF_singleMu_IsoMu24Central
                × inv_N_orig   [× weight_dy_* per il DY]
```

- `inv_N_orig = 1 / Σ genWeight` sull'intero dataset (dai `report_*.json`, prima di ogni taglio).
- `weight_lumi` = 109948.18 pb⁻¹. `weight_xs` = sezione d'urto × BR in pb.
- SF ID/iso: sotto 30 GeV di `mu*_pt` (corretto) si usano Tight ID + Tight iso, sopra Medium ID + Loose iso. Tutti i `weight_mu*_*_Central` sono salvati e l'espressione esatta è nel metadata.
- `weight_gen_sign` = segno di `genWeight`.
- Il MC 2024 è diviso tra le ere 2024–2026 (`OrthogonalEraTag`): circa il 45% degli eventi NanoAOD resta nel 2024, e `weight_lumi` è la lumi del solo 2024.

Rendimento atteso: Σ `weight__Central` sugli eventi selezionati.

## Quale pT entra dove

| passo | pT usato |
|---|---|
| muoni buoni dello skim (pT > 10, \|η\| < 2.4, medium ID, pfIsoId ≥ 2) | **non corretto** (`raw_noCorr`), η NanoAOD |
| trigger matching (HLT_IsoMu24, pT > 26, oggetto id 13, bit 8, ΔR < 0.4) | **non corretto**, direzione non corretta |
| ordinamento mu1/mu2 | corretto (ScaRe + FSR) |
| finestra di massa dello skim 50–200, regioni di massa | corretto |
| tagli del sync 26/20 (`mu_rawpt_sel`, dentro `baseline`) | **non corretto** |
| flag `mu1_pt_g26`, `mu2_pt_g20`, `mu2_pt_g26` | corretto |
| SF ID/iso/trigger: valore | pT **non corretto**, η corretto |
| SF: scelta tight/medium sotto 30 GeV | pT corretto |
| SF di trigger: quale muone | quello con pT corretto più alto, se ha il matching |
| jet selezionati (pT > 25, \|η\| < 4.7, tight ID, fuori veto map e horn, ΔR(j, μ) > 0.4) | pT jet corretto (JEC + JER), ΔR con i muoni nominali |
| b veto (nessun medium, meno di 2 loose tra i jet selezionati) | jet selezionati, \|η\| < 2.5 |
| coppia VBF (pT > 35/25, m_jj > 400, \|Δη\| > 2.5, fuori dall'horn) | pT jet corretto |

Horn 2024: 2.5 ≤ |η| < 3.0 e pT ≤ 50 GeV.

## Metadata (`NAME_metadata.json`)

1. `datasets`: nome, dataset DAS NanoAOD, numero di file skim, voce e valore della sezione d'urto.
2. `normalisation`: lumi integrata e golden JSON; per il MC, l'espressione di `weight__Central` e i valori di `weight_lumi`, `weight_xs`, `inv_N_orig` e Σ genWeight letti dalla tupla.
3. `corrections`: catena di correzioni data/MC e payload (JEC/JER, PU, SF muoni, ScaRe, veto map, b-tag) con i valori dei WP b-tag. È generato da `sync/python/collect_corrections.py`.
4. `selections`: tagli dello skim, filtro della tupla, espressione di ogni flag e lista dei flag salvati.
5. `conventions`: la tabella "Quale pT entra dove", più i jet.

Ci sono anche il commit git, il numero di file modificati, il comando lanciato e la lista dei branch.

Lettura rapida:

```bash
python3 -c "import json; m = json.load(open('NAME_metadata.json')); print(json.dumps(m['normalisation'], indent=1))"
```

## Controlli dopo la produzione

```bash
python3 - $OUT/sync_tuple_loose_VBFHto2Mu_M125_powheg.root <<'EOF'
import sys, ROOT
d = ROOT.RDataFrame("Events", sys.argv[1])
print("eventi", d.Count().GetValue())
print("VBF && Signal_Fit, pesato", d.Filter("VBF && Signal_Fit").Sum("weight__Central").GetValue())
print("pT raw / corretto leading jet", d.Filter("N_SelectedJets > 0")
      .Define("r", "leadingjet_pt_raw / leadingjet_pt").Mean("r").GetValue())
EOF
```

- Il numero di eventi della tupla deve essere uguale all'ultimo step del report skim (`No_extra_electrons`) nel cutflow, perché lo skim ha già 50–200 GeV.
- `add_raw_jet_pt.py` si ferma se un evento non ha un unico evento NanoAOD corrispondente, o se l'allineamento con la tupla non torna. Se arriva in fondo, i branch `*_pt_raw` sono corretti.
- Per i dati, controllare in `NAME_metadata.json` che `datasets` contenga sia Muon0 che Muon1.

## Messaggi noti, innocui

- `Warning in <Snapshot>: A lazy Snapshot action was booked but never triggered`: c'era già con la versione precedente; l'output è completo.
- `Error in <TNetXNGFile::Open>: ... [3001] Required argument not present`, a volte sul MC: intermittente, all'apertura dei file da EOS; non ferma il job.
