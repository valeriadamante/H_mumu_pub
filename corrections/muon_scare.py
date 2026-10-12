#!/usr/bin/env python3
import os
import ROOT
from .general import pog_folder_names,period_names
import correctionlib

correctionlib.register_pyroot_binding()

def apply_muon_scare(
    df,
    config,
    dataset_cfg,
    want_variations,
):
    era = config.get("era")
    period_unc = period_names[era]
    folder_name = pog_folder_names["MUO"][period_unc]

    jsonFile_path = os.path.join(
        os.environ["ANALYSIS_PATH"],
        "corrections",
        "data",
        "MUO",
        "MuonScaRe",
        folder_name,
        "muon_scalesmearing.json"
    )
    jsonFile_path_VXBS = os.path.join(
        os.environ["ANALYSIS_PATH"],
        "corrections",
        "data",
        "MUO",
        "MuonScaRe",
        folder_name,
        "muon_scalesmearing_VXBS.json"
    )

    ## first on raw nano pT

    ROOT.gROOT.ProcessLine(
        f'auto cset = correction::CorrectionSet::from_file("{jsonFile_path}");'
    )
    ROOT.gROOT.ProcessLine(
        f'auto cset_vxbs = correction::CorrectionSet::from_file("{jsonFile_path_VXBS}");'
    )

    headers_dir = os.path.dirname(os.path.abspath(__file__))
    header_path = os.path.join(headers_dir, "MuonScaRe.cc")
    ROOT.gInterpreter.Declare(f'#include "{header_path}"')

    # CENTRAL - Scale (both data and MC)
    df = df.Define(
        'Muon_pt_nano_scale',
        f'''
        ROOT::VecOps::RVec<float> Muon_pt_nano_scale(Muon_pt.size());
        for(size_t Muon_pt_idx = 0 ; Muon_pt_idx < Muon_pt.size(); Muon_pt_idx ++ ){{
            Muon_pt_nano_scale[Muon_pt_idx] = pt_scale(is_data_int, Muon_pt[Muon_pt_idx], Muon_eta[Muon_pt_idx], Muon_phi[Muon_pt_idx], Muon_charge[Muon_pt_idx], false);
        }}
        return Muon_pt_nano_scale;

        '''
    )


    # CENTRAL - Resol (only MC) on top of Scale

    if not dataset_cfg.get("is_data", False):
        df = df.Define(
            'Muon_pt_nano_corr',
            f'''
            ROOT::VecOps::RVec<float> Muon_pt_nano_corr(Muon_pt.size());
            for(size_t Muon_pt_idx = 0 ; Muon_pt_idx < Muon_pt.size(); Muon_pt_idx ++ ){{
                Muon_pt_nano_corr[Muon_pt_idx] = pt_resol(Muon_pt_nano_scale[Muon_pt_idx], Muon_eta[Muon_pt_idx], Muon_phi[Muon_pt_idx], float(Muon_nTrackerLayers[Muon_pt_idx]), event, luminosityBlock, false);
            }}
            return Muon_pt_nano_corr;

            '''
        )
    else:
        df = df.Define('Muon_pt_nano_corr',"Muon_pt_nano_scale")

    # UP/DOWN - Scale, then Resol (only MC)
    if not dataset_cfg.get("is_data", False) and want_variations:
        df = df.Define(
            'Muon_pt_nano_scale_up',
            f'''
            ROOT::VecOps::RVec<float> Muon_pt_nano_scale_up(Muon_pt.size());
            for(size_t Muon_pt_idx = 0 ; Muon_pt_idx < Muon_pt.size(); Muon_pt_idx ++ ){{
                Muon_pt_nano_scale_up[Muon_pt_idx] = pt_scale_var(Muon_pt_nano_corr[Muon_pt_idx], Muon_eta[Muon_pt_idx], Muon_phi[Muon_pt_idx], Muon_charge[Muon_pt_idx], "up", false);
            }}
            return Muon_pt_nano_scale_up;
            '''
        )
        df = df.Define(
            'Muon_pt_nano_scale_down',
            f'''
            ROOT::VecOps::RVec<float> Muon_pt_nano_scale_down(Muon_pt.size());
            for(size_t Muon_pt_idx = 0 ; Muon_pt_idx < Muon_pt.size(); Muon_pt_idx ++ ){{
                Muon_pt_nano_scale_down[Muon_pt_idx] = pt_scale_var(Muon_pt_nano_corr[Muon_pt_idx], Muon_eta[Muon_pt_idx], Muon_phi[Muon_pt_idx], Muon_charge[Muon_pt_idx], "dn", false);
            }}
            return Muon_pt_nano_scale_down;
            '''
        )

        df = df.Define(
                "Muon_pt_nano_res_up",
                f'''
                    ROOT::VecOps::RVec<float> Muon_pt_nano_res_up(Muon_pt.size());
                    for(size_t Muon_pt_idx = 0 ; Muon_pt_idx < Muon_pt.size(); Muon_pt_idx ++ ){{
                        Muon_pt_nano_res_up[Muon_pt_idx] = pt_resol_var(Muon_pt_nano_scale[Muon_pt_idx], Muon_pt_nano_corr[Muon_pt_idx], Muon_eta[Muon_pt_idx], "up", false);
                    }}
                    return Muon_pt_nano_res_up;
                '''
        )

        df = df.Define(
                "Muon_pt_nano_res_down",
                f'''
                    ROOT::VecOps::RVec<float> Muon_pt_nano_res_down(Muon_pt.size());
                    for(size_t Muon_pt_idx = 0 ; Muon_pt_idx < Muon_pt.size(); Muon_pt_idx ++ ){{
                        Muon_pt_nano_res_down[Muon_pt_idx] =  pt_resol_var(Muon_pt_nano_scale[Muon_pt_idx], Muon_pt_nano_corr[Muon_pt_idx], Muon_eta[Muon_pt_idx], "dn", false);
                    }}
                    return Muon_pt_nano_res_down;
                '''
        )

    ## then on bsc pT

    # CENTRAL - Scale (both data and MC)
    df = df.Define(
        'Muon_pt_bsc_scale',
        f'''
        ROOT::VecOps::RVec<float> Muon_pt_bsc_scale(Muon_bsConstrainedPt.size());
        for(size_t Muon_pt_bsc_idx = 0 ; Muon_pt_bsc_idx < Muon_bsConstrainedPt.size(); Muon_pt_bsc_idx ++ ){{
            Muon_pt_bsc_scale[Muon_pt_bsc_idx] = pt_scale(is_data_int, Muon_bsConstrainedPt[Muon_pt_bsc_idx], Muon_eta[Muon_pt_bsc_idx], Muon_phi[Muon_pt_bsc_idx], Muon_charge[Muon_pt_bsc_idx], true);
        }}
        return Muon_pt_bsc_scale;

        '''
    )

    # CENTRAL - Resol (only MC) on top of Scale

    if not dataset_cfg.get("is_data", False):
        df = df.Define(
            'Muon_pt_bsc_corr',
            f'''
            ROOT::VecOps::RVec<float> Muon_pt_bsc_corr(Muon_bsConstrainedPt.size());
            for(size_t Muon_pt_bsc_idx = 0 ; Muon_pt_bsc_idx < Muon_bsConstrainedPt.size(); Muon_pt_bsc_idx ++ ){{
                Muon_pt_bsc_corr[Muon_pt_bsc_idx] = pt_resol(Muon_pt_bsc_scale[Muon_pt_bsc_idx], Muon_eta[Muon_pt_bsc_idx], Muon_phi[Muon_pt_bsc_idx], float(Muon_nTrackerLayers[Muon_pt_bsc_idx]), event, luminosityBlock, true);
            }}
            return Muon_pt_bsc_corr;

            '''
        )
    else:
        df = df.Define('Muon_pt_bsc_corr',"Muon_pt_bsc_scale")

    # UP/DOWN - Scale, then Resol (only MC)
    if not dataset_cfg.get("is_data", False):
        df = df.Define(
            'Muon_pt_bsc_scale_up',
            f'''
            ROOT::VecOps::RVec<float> Muon_pt_bsc_scale_up(Muon_bsConstrainedPt.size());
            for(size_t Muon_pt_bsc_idx = 0 ; Muon_pt_bsc_idx < Muon_bsConstrainedPt.size(); Muon_pt_bsc_idx ++ ){{
                Muon_pt_bsc_scale_up[Muon_pt_bsc_idx] = pt_scale_var(Muon_pt_bsc_corr[Muon_pt_bsc_idx], Muon_eta[Muon_pt_bsc_idx], Muon_phi[Muon_pt_bsc_idx], Muon_charge[Muon_pt_bsc_idx], "up", true);
            }}
            return Muon_pt_bsc_scale_up;
            '''
        )
        df = df.Define(
            'Muon_pt_bsc_scale_down',
            f'''
            ROOT::VecOps::RVec<float> Muon_pt_bsc_scale_down(Muon_bsConstrainedPt.size());
            for(size_t Muon_pt_bsc_idx = 0 ; Muon_pt_bsc_idx < Muon_bsConstrainedPt.size(); Muon_pt_bsc_idx ++ ){{
                Muon_pt_bsc_scale_down[Muon_pt_bsc_idx] = pt_scale_var(Muon_pt_bsc_corr[Muon_pt_bsc_idx], Muon_eta[Muon_pt_bsc_idx], Muon_phi[Muon_pt_bsc_idx], Muon_charge[Muon_pt_bsc_idx], "dn", true);
            }}
            return Muon_pt_bsc_scale_down;
            '''
        )

        df = df.Define(
                "Muon_pt_bsc_res_up",
                f'''
                    ROOT::VecOps::RVec<float> Muon_pt_bsc_res_up(Muon_bsConstrainedPt.size());
                    for(size_t Muon_pt_bsc_idx = 0 ; Muon_pt_bsc_idx < Muon_bsConstrainedPt.size(); Muon_pt_bsc_idx ++ ){{
                        Muon_pt_bsc_res_up[Muon_pt_bsc_idx] = pt_resol_var(Muon_pt_bsc_scale[Muon_pt_bsc_idx], Muon_pt_bsc_corr[Muon_pt_bsc_idx], Muon_eta[Muon_pt_bsc_idx], "up", true);
                    }}
                    return Muon_pt_bsc_res_up;
                '''
        )

        df = df.Define(
                "Muon_pt_bsc_res_down",
                f'''
                    ROOT::VecOps::RVec<float> Muon_pt_bsc_res_down(Muon_bsConstrainedPt.size());
                    for(size_t Muon_pt_bsc_idx = 0 ; Muon_pt_bsc_idx < Muon_bsConstrainedPt.size(); Muon_pt_bsc_idx ++ ){{
                        Muon_pt_bsc_res_down[Muon_pt_bsc_idx] = pt_resol_var(Muon_pt_bsc_scale[Muon_pt_bsc_idx], Muon_pt_bsc_corr[Muon_pt_bsc_idx], Muon_eta[Muon_pt_bsc_idx], "dn", true);
                    }}
                    return Muon_pt_bsc_res_down;
                '''
        )
    return df