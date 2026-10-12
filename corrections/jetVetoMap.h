#pragma once

#include "correction.h"
#include "corrections.h"

namespace correction {
    class JetVetoMapProvider : public CorrectionsBase<JetVetoMapProvider> {
      public:
        JetVetoMapProvider(std::string const& fileName, std::string const& entry_name)
            : corrections_(CorrectionSet::from_file(fileName)) {
            std::cout << "JetVetoMapProvider: init" << std::endl;
            JetVetoMap_value = corrections_->at(entry_name);
        }
        RVecB GetJetVetoMapValues(const RVecLV& Jet_p4) const {
            // Non-zero value for (eta, phi) indicates that the region is vetoed. --> initialization with 0.
            RVecB veto_values(Jet_p4.size(), false);  // Default value is false (NO veto)
            for (int jet_idx = 0; jet_idx < Jet_p4.size(); jet_idx++) {
                float jvm_factor =
                    JetVetoMap_value->evaluate({"jetvetomap", Jet_p4[jet_idx].Eta(), Jet_p4[jet_idx].Phi()});
                veto_values[jet_idx] = static_cast<bool>(jvm_factor);  // Apply the correction factor
            }
            return veto_values;
        }

      private:
        std::unique_ptr<CorrectionSet> corrections_;
        Correction::Ref JetVetoMap_value;
    };

}  // namespace correction

namespace JetIdNewDef{

    RVecB Redefine_Jet_passJetIdTightLepVeto(const RVecLV& Jet_p4, const RVecB& Jet_passJetIdTight, const RVecF& Jet_muEF, const RVecF& Jet_chEmEF){
        RVecB Jet_passJetIdTightLepVeto(Jet_p4.size(),false);
        for(size_t jet_index = 0; jet_index < Jet_p4.size(); jet_index ++){
            if(abs(Jet_p4[jet_index].eta()) <= 2.7){
                Jet_passJetIdTightLepVeto[jet_index] = Jet_passJetIdTight[jet_index] && (Jet_muEF[jet_index] < 0.8) && (Jet_chEmEF[jet_index] < 0.8);
            }
            else{
                Jet_passJetIdTightLepVeto[jet_index] = Jet_passJetIdTight[jet_index];
            }
        }
        return Jet_passJetIdTightLepVeto;
    }
    RVecB RedefineJet_passJetIdTight_v12(const RVecLV& Jet_p4, const RVecF& Jet_neHEF, const RVecF& Jet_neEmEF, const RVecUC& Jet_jetId){
        RVecB Jet_passJetIdTight(Jet_p4.size(),false);
        for(size_t jet_index = 0; jet_index < Jet_p4.size(); jet_index ++){
            if(abs(Jet_p4[jet_index].eta()) <= 2.7){
                Jet_passJetIdTight[jet_index] = Jet_jetId[jet_index] & (1 << 1);
            }
            else if(abs(Jet_p4[jet_index].eta()) > 2.7 && abs(Jet_p4[jet_index].eta()) <= 3.){
                Jet_passJetIdTight[jet_index] = Jet_jetId[jet_index] & (1 << 1) && (Jet_neHEF[jet_index] < 0.99);
            }
            else if(abs(Jet_p4[jet_index].eta()) > 3.){
                Jet_passJetIdTight[jet_index] = Jet_jetId[jet_index] & (1 << 1) && (Jet_neEmEF[jet_index] < 0.99);
            }
        }
        return Jet_passJetIdTight ;
    }

    RVecB RedefineJet_passJetIdTight_v13(const RVecLV& Jet_p4, const RVecF& Jet_neHEF, const RVecF& Jet_neEmEF, const RVecF& Jet_chHEF, const RVecF& Jet_chMultiplicity, const RVecF& Jet_neMultiplicity ){
        RVecB Jet_passJetIdTight(Jet_p4.size(),false);
        for(size_t jet_index = 0; jet_index < Jet_p4.size(); jet_index ++){
            if(abs(Jet_p4[jet_index].eta()) <= 2.6){
                Jet_passJetIdTight[jet_index] = (Jet_neHEF[jet_index] < 0.99) && (Jet_neEmEF[jet_index] < 0.9) && (Jet_chMultiplicity[jet_index]+Jet_neMultiplicity[jet_index] > 1) && (Jet_chHEF[jet_index] > 0.01) && (Jet_chMultiplicity[jet_index] > 0);
            }
            else if(abs(Jet_p4[jet_index].eta()) > 2.6 && abs(Jet_p4[jet_index].eta()) <= 2.7){
                Jet_passJetIdTight[jet_index] = (Jet_neHEF[jet_index] < 0.90) && (Jet_neEmEF[jet_index] < 0.99);
            }
            else if(abs(Jet_p4[jet_index].eta()) > 2.7 && abs(Jet_p4[jet_index].eta()) <= 3.){
                Jet_passJetIdTight[jet_index] = (Jet_neHEF[jet_index] < 0.99);
            }
            else if(abs(Jet_p4[jet_index].eta()) > 3.){
                Jet_passJetIdTight[jet_index] = (Jet_neMultiplicity[jet_index] >= 2) && (Jet_neEmEF[jet_index] < 0.4);
            }
        }
        return Jet_passJetIdTight ;
    }
}