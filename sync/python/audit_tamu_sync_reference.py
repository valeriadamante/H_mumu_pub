#!/usr/bin/env python3
"""Audit the new TAMU Run2024H loose tuple against the 2026-09-25 Pisa-sync tuple."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import uproot

FILES={
 'reference_filippo_20260925':Path('/eos/user/v/vdamante/H_mumu/Sync/sync_tuple_Run2024H_filippo_20260925.root'),
 'new_TAMU':Path('/eos/cms/store/group/phys_higgs/cmshmm/run3_sync/TAMU/TAMU-Run2024H.root')}
FLAGS=['HLT_IsoMu24','Event_HasTriggerMatching_singleMu','muons_OS','base_sel_muons','electron_veto','jets_btag_clean_sel','baseline','VBF','ggF','Z_sideband','H_sideband','Signal_Fit','Signal_ext','mass_inclusive']
VALUES=['m_mumu_FSR_corr','mu1_pt_FSR_corr','mu2_pt_FSR_corr','mu1_eta','mu2_eta','N_SelectedJets','leadingjet_pt','subleadingjet_pt','dijet_ls_mass']

def t(path): return uproot.open(path)['Events']
def keys(tree):
 a=tree.arrays(['run','luminosityBlock','event'],library='np')
 k=np.rec.fromarrays([a['run'].astype('u8'),a['luminosityBlock'].astype('u8'),a['event'].astype('u8')],names='run,lumi,event')
 return np.unique(k,return_index=True)
def stats(tree):
 present=[x for x in FLAGS if x in tree.keys()]; out={}
 for a in tree.iterate(present,step_size='200 MB',library='np'):
  for x in present: out[x]=out.get(x,0)+int(a[x].astype(bool).sum())
 return out
def matched_values(tree,branches,indices):
 a=tree.arrays(branches,library='np'); return {x:np.asarray(a[x])[indices] for x in branches}
def main():
 p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True); a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
 trees={n:t(f) for n,f in FILES.items()}; ku={}; first={}; summary={'files':{n:str(f) for n,f in FILES.items()},'entries':{},'unique':{},'duplicates':{},'flags':{}}
 for n,tr in trees.items():
  k,i=keys(tr); ku[n]=k; first[n]=i; summary['entries'][n]=tr.num_entries; summary['unique'][n]=len(k); summary['duplicates'][n]=tr.num_entries-len(k); summary['flags'][n]=stats(tr)
 old,new=ku; common=np.intersect1d(ku[old],ku[new]); io=np.searchsorted(ku[old],common); inn=np.searchsorted(ku[new],common); summary['common']=len(common)
 summary['only_reference']=len(ku[old])-len(common); summary['only_new']=len(ku[new])-len(common)
 bo=[x for x in [*FLAGS,*VALUES] if x in trees[old].keys() and x in trees[new].keys()]
 vo=matched_values(trees[old],bo,first[old][io]); vn=matched_values(trees[new],bo,first[new][inn]); summary['matched']={}
 for x in bo:
  xo,xn=vo[x],vn[x]
  if xo.dtype.kind in 'bui' and xn.dtype.kind in 'bui':
   summary['matched'][x]={'disagree':int((xo!=xn).sum()),'fraction_disagree':float((xo!=xn).mean())}
  else:
   d=xn.astype(float)-xo.astype(float); finite=np.isfinite(d)
   summary['matched'][x]={'different_gt_1e-6':int((np.abs(d[finite])>1e-6).sum()),'fraction_different':float((np.abs(d[finite])>1e-6).mean()),'median_new_minus_reference':float(np.median(d[finite])),'p99_abs_difference':float(np.quantile(np.abs(d[finite]),.99))}
 (a.output/'audit.json').write_text(json.dumps(summary,indent=2)+'\n')
 lines=['# TAMU Run2024H tuple audit','', 'Comparison of the new shared TAMU tuple with the established `filippo_20260925` reference tuple.','', '## Event coverage','', '| tuple | entries | unique | duplicates |','|---|---:|---:|---:|']
 for n in FILES: lines.append(f'| {n} | {summary["entries"][n]:,} | {summary["unique"][n]:,} | {summary["duplicates"][n]:,} |')
 lines += ['',f'- common: **{len(common):,}**',f'- reference only: **{summary["only_reference"]:,}**',f'- new TAMU only: **{summary["only_new"]:,}**','', '## Selection and category counts','', '| flag | reference | new TAMU | new/reference |','|---|---:|---:|---:|']
 for flag in FLAGS:
  x=summary['flags'][old].get(flag); y=summary['flags'][new].get(flag)
  if x is not None and y is not None: lines.append(f'| `{flag}` | {x:,} | {y:,} | {y/x:.4f} |' if x else f'| `{flag}` | 0 | {y:,} | n/a |')
 lines += ['', '## Agreement on common events','', '| branch | diagnostic | value |','|---|---|---:|']
 for branch,v in summary['matched'].items():
  if 'fraction_disagree' in v: lines.append(f'| `{branch}` | fraction with different boolean/integer value | {v["fraction_disagree"]:.4%} |')
  else: lines.append(f'| `{branch}` | fraction differing by >1e-6 | {v["fraction_different"]:.4%} |')
 lines += ['', '## Interpretation guide','', '- Large event-set differences with stable values on common events indicate input coverage or tuple preselection, not reconstruction.', '- Differences in muon/dimuon values on common events point to calibration/FSR or pairing changes.', '- Stable muons with changing `N_SelectedJets`, jet values, `VBF`, or `ggF` point to JEC/JER, cleaning, veto maps, ordering, or category definitions.', '- The JSON contains median and 99th-percentile numerical differences for follow-up.','']
 (a.output/'SUMMARY.md').write_text('\n'.join(lines))
if __name__=='__main__': main()
