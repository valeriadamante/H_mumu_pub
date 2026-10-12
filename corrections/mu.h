#pragma once

#include <string>
#include <iostream>

namespace correction {

    // A simple, clean stateless evaluator function
    inline double getMuonSF_simple(
        const std::string& cset_key,
        const std::string& systematic, // "nominal", "systup", or "systdown"
        float pt,
        float eta,
        float pfRelIso04_all,
        bool tightId,
        float tkRelIso,
        // bool highPtId,
        bool mediumId,
        bool looseId,
        bool trg_matching,
        bool trg_path,
        int event
    ) {

        // 3. Keep your custom object selection gating if needed.
        // This ensures you only apply corrections where the object passes the selection.
        const bool LooseIso = (pfRelIso04_all < 0.25);
        const bool MediumIso = (pfRelIso04_all < 0.2);
        const bool TightIso = (pfRelIso04_all < 0.15);
        const bool LooseRelTrkIso = (tkRelIso < 0.25);
        const bool MediumRelTrkIso = (tkRelIso < 0.2);
        const bool TightRelTrkIso = (tkRelIso < 0.15);

        // ======================
        // RECO
        // ======================
        if (cset_key=="NUM_TrackerMuons_DEN_genTracks")
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_GlobalMuons_DEN_genTracks")
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        // ======================
        // ID (genTracks / TrackerMuons)
        // ======================
        if (cset_key=="NUM_LooseID_DEN_genTracks" && looseId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LooseID_DEN_TrackerMuons" && looseId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_MediumID_DEN_genTracks" && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_MediumID_DEN_TrackerMuons" && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_MediumPromptID_DEN_genTracks" && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_MediumPromptID_DEN_TrackerMuons" && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_TightID_DEN_genTracks" && tightId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_TightID_DEN_TrackerMuons" && tightId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_SoftID_DEN_genTracks")
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_SoftID_DEN_TrackerMuons")
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        // ======================
        // ISO – PF
        // ======================
        if (cset_key=="NUM_LoosePFIso_DEN_LooseID" && LooseIso && looseId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LoosePFIso_DEN_MediumID" && LooseIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LoosePFIso_DEN_MediumPromptID" && LooseIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LoosePFIso_DEN_TightID" && LooseIso && tightId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_TightPFIso_DEN_MediumID" && TightIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_TightPFIso_DEN_MediumPromptID" && TightIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_TightPFIso_DEN_TightID" && TightIso && tightId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        // ======================
        // ISO – RelIso (PF + tk)
        // ======================

        if (cset_key=="NUM_LooseRelIso_DEN_LooseID" && LooseIso && looseId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LooseRelIso_DEN_MediumID" && LooseIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LooseRelIso_DEN_MediumPromptID" && LooseIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_LooseRelIso_DEN_TightIDandIPCut" && LooseIso && tightId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key=="NUM_TightRelIso_DEN_MediumID" && TightIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_TightRelIso_DEN_MediumPromptID" && TightIso && mediumId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key=="NUM_TightRelIso_DEN_TightIDandIPCut" && TightIso && tightId)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});

        if (cset_key =="NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium" && (mediumId || looseId) && (LooseIso || MediumIso) && trg_matching && trg_path)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        if (cset_key =="NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight" && tightId && TightIso && trg_matching && trg_path)
            return cset->at(cset_key)->evaluate({eta, pt,systematic});
        return 1.0;


    }
}