# Sync MC con Pisa: VBFH → μμ, 2024, regione di segnale VBF

Data: 2026-10-07. Campione: `VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8`
(RunIII2024Summer24NanoAODv15, v2 + ext1). Nostro nome: `VBFHto2Mu_M125_powheg`.

## Input

| | TAMU | Pisa (Filippo) |
|---|---|---|
| file | 190 NanoAOD (v2 + ext1) | gli stessi 190 (`Runs` con 190 voci) |
| eventi iniziali | 11 190 000 | 11 190 000 |
| ntupla | `/eos/user/v/vdamante/H_mumu/Sync/sync_tuple_VBFHto2Mu_M125_powheg_filippo_20260925.root` | `/eos/home-f/fcattafe/Hmm_sync_valeria/hmm_2024Summer24_VBF-Hto2Mu_M-125_powheg_20261005T124327_VBF_SignalRegion_diagnostic_ntuples.root` |
| selezione | skim speciale `filippo_20260925` (definizioni Pisa, `config/Run3_2024/selections_sync_filippo_20260925.yaml`), `VBF && Signal_Fit` | tutti i passi fino a `HiggsMass` |
| finestra di massa | 115 < m_μμ < 135 GeV | `HiggsMass`: 115–135 GeV (verificato sulla ntupla) |
| xsec | 4.078 × 2.176e-4 = 8.8737e-4 pb | 8.8819e-4 pb |
| lumi 2024 | 109 948.18 pb⁻¹ | 109 341 pb⁻¹ |
| Σ genEventSumw | 46 762 888 (stessi file) | 46 762 888 |

Il suo ultimo passo, `MC_selection`, tiene in modo casuale circa il 50% degli eventi.
Non è correlato con kFold, con la parità di `event` né con la lumisection, e non ha
un fattore di peso. Non si può riprodurre, quindi il confronto si ferma a `HiggsMass`.

I suoi tagli in pT usano il pT NanoAOD (`LeadingMuon_pt`); `*_ptToUse` è il pT
corretto. Il minimo di `LeadingMuon_ptToUse` sugli eventi selezionati è 25.4 GeV.

## Confronto evento per evento (senza pesi)

| | eventi |
|---|---|
| TAMU | 2 294 633 |
| Pisa | 2 284 560 |
| in comune | 2 257 365 (98.4% dei nostri, 98.8% dei suoi) |
| solo TAMU | 37 268 |
| solo Pisa | 27 195 |

Le chiavi run:lumi:event duplicate perché riusate da v2 ed ext1 sono 69 da noi e 68 da
lui; nel conteggio sopra sono contate una volta sola.

Solo TAMU: il suo primo passo fallito.

| passo | eventi |
|---|---|
| TwoJets | 23 714 |
| DijetMass | 7 841 |
| BJetVeto | 3 742 |
| LeadingJetPt | 1 035 |
| HiggsMass | 757 |
| TwoMuons | 165 |
| DijetDeltaEta | 9 |
| AtLeastOneTrgMuon | 5 |

Solo Pisa: il nostro primo passo fallito.

| passo | eventi |
|---|---|
| pT dei due getti leading (35/25 GeV) | 11 115 |
| m_jj > 400 GeV | 8 722 |
| b-veto | 6 297 |
| veto elettroni | 534 |
| finestra di massa | 478 |
| non nella nostra tupla | 28 |
| Δη_jj > 2.5 | 21 |

- **Getti.** È la differenza principale: noi abbiamo l'1.2% in più di eventi con ≥ 2 getti
  (3 990 072 contro 3 943 668). Gli eventi migrano in entrambe le direzioni attorno a
  25/35 GeV e a m_jj = 400 GeV. È un'ipotesi da verificare: lo smearing JER in MC è
  casuale, quindi seed diversi danno pT diversi vicino alle soglie. Da chiedere a
  Filippo: seed JER, jet ID, PU ID.
- **Veto elettroni.** Il suo è diverso dal nostro: dei 534 eventi persi, 490 si salvano
  chiedendo pT > 20 GeV e 288 usando il WP80.
- **Massa.** Il rapporto noi/Pisa di m_μμ sale da 0.99 a 115 GeV fino a 1.04 a 135 GeV,
  coerente con i 757 + 478 eventi che attraversano i bordi della finestra. È una piccola
  differenza nella scala dei muoni corretti o nell'FSR; il pT corretto del leading muon
  è entro lo 0.5%.

Cutflow nell'ordine di Pisa: sue entries (dalla tabella) e nostri conteggi del passo
corrispondente (`*_cutflow_vbf_Z_CR_H_SB.csv`). Fino al trigger matching l'ordine
degli step è diverso, quindi lì il confronto è solo indicativo.

| Pisa | entries Pisa | TAMU | conteggio TAMU |
|---|---|---|---|
| NoiseFilters | 11 180 406 | MET_filters | 11 180 406 |
| isTriggered + AtLeastOneTrgMuon | 9 976 479 | Trigger_matching_for_singleMu | 10 078 151 |
| LeadingMuonCut (≈ due muoni, OS, 26/20) | 7 178 647 | mu2_pt_g20 | 7 172 547 |
| TwoJets | 3 943 668 | jets_ge2 | 3 990 072 |
| LeadingJetPt | 3 871 120 | filippo_LeadingJetPt | 3 913 826 |
| DijetMass | 2 518 049 | filippo_DijetMass | 2 530 115 |
| DijetDeltaEta | 2 479 967 | filippo_DijetDeltaEta | 2 491 998 |
| BJetVeto | 2 367 333 | jets_btag_clean_sel | 2 376 355 |
| EleVeto | 2 365 959 | electron_veto | 2 374 304 |
| AdditionalMuonVeto | 2 365 133 | (nello skim) | 2 374 304 |
| HiggsMass | 2 284 628 | Signal_Fit | 2 294 702 |

## Pesi

### I pesi di Pisa

Ricavati dalla ntupla: il prodotto dei `*_weightFactor` riproduce esattamente
`HiggsMass_cumulativeWeight`. I `cumulativeWeight` non sono divisi per la somma dei pesi.

| suo passo | fattore sugli eventi selezionati | significato (da confermare con Filippo) |
|---|---|---|
| BASESELECTION | ±4.184 | genWeight |
| GoldenJson | 0.27–73, media 1.00 | pileup |
| TwoMuons | 0.89–1.05 | SF ID × isolamento dei due muoni |
| LeadingMuonCut | 0.94–1.03 | SF di trigger |
| TwoJets, BJetVeto | sempre 1 | nessun SF su getti o b-tag |

### I nostri pesi

Calcolati da `sync/python/sync_mc_weights.py` con le funzioni dell'analisi, cioè la
voce `Central` di `config/Run3_2024/systematics.yaml` senza xsec, lumi e 1/Σw:

```
genWeight × weight_pu
  × Π_i [(mu_i_pt < 30) ? TightPFIso_TightID × TightID_Trk : LoosePFIso_MediumID × MediumID_Trk]
  × weight_TrgSF_singleMu_IsoMu24 (SF IsoMu24 Medium del muone leading, se ha il trigger match)
```

- **PU.** Tabella di `corrections/pu.py` (`Collisions24_BCDEFGHI_goldenJSON`, valutata a
  floor(`Pileup_nTrueInt`)). `Pileup_nTrueInt` è letto dai 190 file NanoAOD e associato
  per run:lumi:event; 510 eventi hanno chiavi riusate tra v2 ed ext1.
- **SF dei muoni.** `corrections/mu.py` e `muon_Z.json.gz`, cartella `latest` di
  `Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15`. Sono valutati a
  `mu_pt_raw_noCorr`, con η e pT del muone = bsc + ScaRe + FSR, come nello skim di
  produzione.
- **Trigger.** `common/trigger_weights.py`.

Il nostro peso è quindi la stessa cosa del suo: genWeight × PU × SF ID/iso × SF trigger.
Gli altri pesi dell'analisi (correzioni DY) non riguardano questo campione.

### Confronto fattore per fattore sugli eventi in comune

Il confronto usa 2 257 294 eventi; quelli con chiavi duplicate sono esclusi.

| fattore | media TAMU | media Pisa | Σ TAMU / Σ Pisa | rapporto evento per evento p1 / p50 / p99 | uguali entro 0.1% |
|---|---|---|---|---|---|
| genWeight | 4.1832 | 4.1833 | 1.0000 | 1.0000 / 1.0000 / 1.0000 | 100% |
| pileup | 0.9997 | 0.9997 | 1.0000 | 1.0000 / 1.0000 / 1.0000 | 100% |
| ID × iso (2 muoni) | 0.9849 | 0.9950 | 0.9899 | 0.922 / 0.992 / 1.066 | 3.2% |
| trigger | 0.9783 | 0.9872 | 0.9910 | 0.917 / 0.983 / 1.058 | 1.3% |
| totale | 4.0293 | 4.1077 | 0.9809 | 0.911 / 0.975 / 1.030 | 0.9% |

- **genWeight e PU.** Identici evento per evento.
- **SF dei muoni.** Sono diversi. Nessuna combinazione riproduce i suoi `TwoMuons` e
  `LeadingMuonCut`: ho provato le chiavi Medium/Loose, Tight e Medium-ID only, i pT
  NanoAOD, bsc, ScaRe e FSR, η con segno o |η|, tutte le versioni del `muon_Z.json.gz`
  2024 su cvmfs (2025-08-14 → latest) e anche i file 2022 e 2023. La differenza
  dipende da η: con entrambi i muoni con |η| < 0.9 il rapporto mediano noi/Pisa
  dell'ID × iso è 0.994, con entrambi tra 2.1 e 2.4 è 0.965. Nei plot pesati, η(μ1)
  va da circa −3% nel barrel a +2% negli endcap. Da chiedere a Filippo: quale file e
  quali chiavi usa per ID, isolamento e trigger, e a quale pT li valuta.
- **pT < 30 GeV.** La nostra analisi usa le chiavi Tight sotto i 30 GeV; con Medium/Loose
  anche lì il disaccordo resta (rapporto mediano 0.992).

### Yield nella regione VBF, 115–135 GeV

Normalizzazione: Σw / Σ genEventSumw × xsec × lumi, prima di `MC_selection`.
Dove non indicato, xsec e lumi sono quelle di Pisa.

| | eventi |
|---|---|
| TAMU, solo genWeight | 19.94 |
| TAMU, pesi completi | 19.23 |
| Pisa, solo genWeight | 19.85 |
| Pisa, `HiggsMass_cumulativeWeight` | 19.51 |
| TAMU, pesi completi, nostra xsec e lumi | 19.32 |

TAMU/Pisa vale +0.4% con il solo genWeight (selezione) e −1.4% con i pesi completi;
la differenza tra i due viene dagli SF dei muoni.

**Correzione rispetto a ieri.** Avevo segnalato un fattore 2 nella sua tabella; non è
un errore. La sua riga finale (`MC_selection`, 0.00017869 pb × 109 341 = 19.54 eventi)
coincide con la somma di `HiggsMass_cumulativeWeight` normalizzata alla sumw totale
(19.51). La sua tabella normalizza per la sumw della metà tenuta da `MC_selection`:
per questo al passo 0 vale 2 × xsec e il numero finale è giusto. I 9.7 eventi che
avevo calcolato erano la somma dei pesi dopo `MC_selection` divisa per la sumw totale,
cioè la metà del valore corretto.

## Plot

Sono 17 distribuzioni TAMU contro Pisa, con il rapporto; l'incertezza del rapporto
conta solo gli eventi non in comune.

- solo genWeight: `/eos/user/v/vdamante/www/H_mumu/skim_v4/studies/sync/filippo_mc_vbfh/`
- pesi completi (`weight_analysis` contro `HiggsMass_cumulativeWeight`):
  `/eos/user/v/vdamante/www/H_mumu/skim_v4/studies/sync/filippo_mc_vbfh_allweights/`

## File prodotti

In `/eos/user/v/vdamante/H_mumu/Sync/confronto_filippo_mc_vbfh/`:

- `riassunto.csv`, `solo_nostri.csv`, `solo_filippo.csv`: confronto evento per evento,
  con il primo passo fallito dall'altra parte.
- `weights_VBF_Signal_Fit.root`: i nostri pesi per evento (`genWeight`, `weight_pu`,
  `weight_muIDIso`, `weight_trigger`, `weight_analysis`, `Pileup_nTrueInt`).
- `pileup_nanoaod_VBFHto2Mu_M125_powheg.npz`: `Pileup_nTrueInt` dai NanoAOD.
- `weights.md`, `weight_ratios.csv`: il confronto dei pesi.

Accanto alla tupla ci sono i cutflow CSV: `sync_tuple_VBFHto2Mu_M125_powheg_filippo_20260925_*.csv`.

## Comandi

```bash
source env.sh
FIL='/eos/home-f/fcattafe/Hmm_sync_valeria/hmm_2024Summer24_VBF-Hto2Mu_M-125_powheg_*_VBF_SignalRegion_diagnostic_ntuples.root'
OUT=$HMUMU_SYNC/confronto_filippo_mc_vbfh
python3 sync/python/special_sync_submit.py --variant filippo_20260925 --dataset VBFHto2Mu_M125_powheg
python3 sync/python/special_sync_merge.py --tag filippo_20260925 --dataset VBFHto2Mu_M125_powheg
python3 sync/python/compare_filippo_diagnostic.py --dataset VBFHto2Mu_M125_powheg --region Signal_Fit \
    --filippo-step HiggsMass --filippo-input "$FIL" --output-dir $OUT
python3 sync/python/sync_mc_weights.py --dataset VBFHto2Mu_M125_powheg --output $OUT/weights_VBF_Signal_Fit.root
python3 sync/python/compare_filippo_weights.py --ours-weights $OUT/weights_VBF_Signal_Fit.root \
    --filippo-input "$FIL" --output-dir $OUT
python3 sync/python/plot_filippo_comparison.py --dataset VBFHto2Mu_M125_powheg --region Signal_Fit \
    --filippo-step HiggsMass --filippo-input "$FIL" --weight analysis \
    --ours-weights $OUT/weights_VBF_Signal_Fit.root --filippo-weight HiggsMass_cumulativeWeight \
    --output-dir $HMUMU_WWW/skim_v4/studies/sync/filippo_mc_vbfh_allweights
```

## Domande per Filippo

1. Quale `muon_Z.json` (versione/data) e quali chiavi usa per ID, isolamento e trigger;
   a quale pT e η li valuta; se il trigger SF è sul solo muone leading.
2. Seed dello smearing JER, jet ID e PU ID: sono la causa probabile delle migrazioni
   attorno alle soglie dei getti.
3. Definizione del veto elettroni: pT > 20 GeV? WP80?
4. Correzione della scala dei muoni e FSR usate per `Dimuon_mass`.
5. `MC_selection`: è uno split train/test? E quale lumi usa (109 341 contro i nostri
   109 948 pb⁻¹)?
