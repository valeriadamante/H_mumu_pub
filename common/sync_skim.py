"""Nominal, unweighted event-level synchronization outputs."""
import json
from pathlib import Path

from common.add_vars import (
    DefineSelections, GetAllMuonsObservablesNew, SelectedJetObservablesDef,
    SoftJetCollectionCleaningInVBF, VBFJetMuonsObservablesDef, VBFJetObservablesDef,
)


def regions_for_sample(is_data):
    return ['H_sideband'] if is_data else ['H_sideband', 'Signal_Fit']


def json_value(value):
    if hasattr(value, 'tolist'):
        return json_value(value.tolist())
    if hasattr(value, 'item'):
        return value.item()
    if not isinstance(value, (str, bytes)) and hasattr(value, '__iter__'):
        return [json_value(item) for item in value]
    return value


def prepare_sync(df, selections, is_data, category):
    # Same order as common.prepare_rdf: the sync must report exactly the
    # high-level quantities used by the analysis, not a redefinition of them.
    df = SelectedJetObservablesDef(df)
    df = VBFJetObservablesDef(df)
    df = GetAllMuonsObservablesNew(df)
    df = VBFJetMuonsObservablesDef(df)
    df = SoftJetCollectionCleaningInVBF(df)
    df = DefineSelections(df, selections)
    df = df.Define('Sync_Jet_pt_raw', 'Jet_pt_nocorr * (1.f - Jet_rawFactor)')
    df = df.Define('Sync_SelectedJet_pt_raw', 'Take(Sync_Jet_pt_raw, SelectedJet_idx)')
    df = df.Filter(category, 'Sync analysis selection: ' + category)
    regions = regions_for_sample(is_data)
    df = df.Filter(' || '.join(regions), 'Sync mass regions: ' + ', '.join(regions))
    extra = [category, *regions, 'Sync_Jet_pt_raw', 'Sync_SelectedJet_pt_raw', 'Jet_rawFactor']
    return df, extra


def export_sync(df, selections, is_data, category, output_dir, metadata):
    import ROOT
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    available = {str(c) for c in df.GetColumnNames()}
    columns = ['run', 'luminosityBlock', 'event', 'm_mumu', 'N_SelectedJets',
               'muons_OS', 'mu_pt_sel', 'mu_pt_trg_sel', category,
               'Jet_pt', 'Jet_pt_nocorr', 'Jet_rawFactor', 'Sync_Jet_pt_raw',
               'Jet_eta', 'Jet_phi', 'Jet_mass', 'Jet_jetId', 'Jet_vetoMap',
               'Jet_IsInsideHorn', 'Jet_preSel', 'goodJet', 'SelectedJet_idx',
               'SelectedJet_pt', 'SelectedJet_eta', 'SelectedJet_phi', 'SelectedJet_mass',
               'Sync_SelectedJet_pt_raw', 'SelectedJetTagSel', 'HasVBF',
               'VBFJetIdx_1', 'VBFJetIdx_2', 'HLT_IsoMu24']
    for leg in (1, 2):
        columns += [f'mu{leg}_{field}' for field in (
            'idx', 'pt', 'eta', 'phi', 'mass', 'charge', 'pt_raw_noCorr',
            'pt_raw_corr', 'pt_FSR_corr', 'mediumId', 'tightId', 'pfIsoId',
            'HasTriggerMatching_singleMu', 'GenMatched', 'genPartIdx', 'genPartFlav')]
    # Dimuon variables (BSC+KIT corrected muons and recovered FSR) and
    # high-level inputs of the ggH/VBF training, so the sync covers the same
    # quantities used by the DNN and not only the objects.
    columns += ['pt_mumu', 'eta_mumu', 'phi_mumu', 'y_mumu', 'dR_mumu',
                'cosTheta_CS', 'phi_CS', 'R_pt', 'minDeltaEtaSigned', 'minDeltaPhi',
                'Zeppenfeld_Var', 'pt_centrality', 'm_jj', 'm_jj_ls',
                'delta_eta_jj', 'delta_eta_jj_ls', 'pt_vbfj1j2', 'era_code',
                'SoftActivityJetHT', 'nSoftActivityJet',
                'SoftJetCleanedActivity_N', 'SoftJetCleanedActivity_ptSum',
                'SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_N',
                'SoftJetActivity_NoOverlapWithMuonsAndEtaCleaning_ptSum']
    for leg in (1, 2):
        columns += [f'vbfjet{leg}_{field}' for field in
                    ('pt', 'eta', 'phi', 'mass', 'y', 'btagPNetQvG')]
    # Weights for the MC synchronization; absent in data and therefore filtered out.
    if not is_data:
        columns += ['genWeight', 'weight_Central', 'weight__Central',
                    'weight_pu_Central', 'puWeight']
        columns += [c for c in sorted(available)
                    if c.startswith(('weight_mu1_', 'weight_mu2_')) and c.endswith('_Central')]
    columns = list(dict.fromkeys(c for c in columns if c in available))
    summary = dict(metadata, category=category, counts={}, regions={},
                   distribution_weights='none (event counts)', event_key=['run', 'luminosityBlock', 'event'],
                   event_columns=columns)
    hist_file = ROOT.TFile.Open(str(out / 'distributions.root'), 'RECREATE')
    if not hist_file or hist_file.IsZombie():
        raise RuntimeError('Cannot create synchronization histograms')
    plots = [('m_mumu', 100, 50., 200.), ('mu1_pt', 100, 0., 200.),
             ('mu2_pt', 100, 0., 200.), ('mu1_eta', 48, -2.4, 2.4),
             ('mu2_eta', 48, -2.4, 2.4), ('N_SelectedJets', 15, -0.5, 14.5),
             ('SelectedJet_pt', 100, 0., 300.), ('SelectedJet_eta', 94, -4.7, 4.7)]
    # Book every action before execution so ROOT reads the NanoAOD only once.
    actions = {}
    for region in regions_for_sample(is_data):
        region_df = df.Filter(region)
        arrays = region_df.AsNumpy(columns, lazy=True)
        histograms = [(name, region_df.Histo1D((name, name, bins, low, high), name))
                      for name, bins, low, high in plots if name in available]
        actions[region] = (arrays, histograms)
    for region, (array_result, histograms) in actions.items():
        arrays = array_result.GetValue()
        order = sorted(range(len(arrays['event'])), key=lambda i: tuple(
            int(arrays[c][i]) for c in ('run', 'luminosityBlock', 'event')))
        ids = []
        with (out / f'{region}_events.jsonl').open('w') as details:
            for i in order:
                row = {c: json_value(arrays[c][i]) for c in columns}
                details.write(json.dumps(row, allow_nan=False) + '\n')
                ids.append(f"{row['run']}:{row['luminosityBlock']}:{row['event']}")
        (out / f'{region}_event_ids.txt').write_text(''.join(key + '\n' for key in ids))
        summary['counts'][region] = len(order)
        summary['regions'][region] = {
            'expression': selections['masses_regions'][region]['expression'].format(mu_suff='', tot_suff='', jet_suff=''),
            'duplicate_event_ids': len(ids) - len(set(ids)),
        }
        directory = hist_file.mkdir(region)
        for name, hist in histograms:
            directory.cd()
            hist.GetValue().Write()
    hist_file.Close()
    summary['selected_events'] = sum(summary['counts'].values())
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('[SYNC] Selected event counts:', summary['counts'])
