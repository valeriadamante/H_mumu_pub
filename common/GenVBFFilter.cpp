#pragma once

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <vector>

using RVecB = ROOT::VecOps::RVec<Bool_t>;
using RVecI = ROOT::VecOps::RVec<Int_t>;
using RVecUS = ROOT::VecOps::RVec<UShort_t>;
using RVecF = ROOT::VecOps::RVec<Float_t>;

bool genVBFFilter_func(const RVecF &GenJet_pt, const RVecF &GenJet_eta, const RVecF &GenJet_phi, const RVecF &GenJet_mass, const RVecF &GenPart_pt, const RVecF &GenPart_eta, const RVecF &GenPart_phi, const RVecF &GenPart_mass, const RVecI &GenPart_pdgId, const RVecUS &GenPart_statusFlags)
{
    std::vector<ROOT::Math::PtEtaPhiMVector> GenLepton;
    for (size_t i = 0; i < GenPart_pt.size(); ++i)
    {
        if (std::abs(GenPart_pdgId[i]) == 11 || std::abs(GenPart_pdgId[i]) == 13 || std::abs(GenPart_pdgId[i]) == 15)
        { // electron or muon or tau
            if ((GenPart_statusFlags[i] & (1 << 7)) != 0)
            { // check isHardProcess bit
                GenLepton.emplace_back(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
            }
        }
    }

    const float minPt = 0.0;
    const float minEta = -99999.0;
    const float maxEta = 99999.0;
    const float minLeadingJetsInvMass = 300.0;
    const float deltaRNoLep = 0.3;

    std::vector<ROOT::Math::PtEtaPhiMVector> GenJetsWithoutLeptons;
    for (size_t i = 0; i < GenJet_pt.size(); ++i)
    {
        if (GenJet_pt[i] > minPt && GenJet_eta[i] > minEta && GenJet_eta[i] < maxEta)
        {
            bool jetWhitoutLep = true;
            ROOT::Math::PtEtaPhiMVector jet(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
            for (const auto &lepton : GenLepton)
            {
                if (ROOT::Math::VectorUtil::DeltaR(jet, lepton) < deltaRNoLep)
                {
                    jetWhitoutLep = false;
                    break;
                }
            }
            if (jetWhitoutLep)
            {
                GenJetsWithoutLeptons.push_back(jet);
            }
        }
    }

    // std::cout << "Number of GenJets: " << GenJet_pt.size() << std::endl;
    // std::cout << "Number of GenJetsWithoutLeptons: " << GenJetsWithoutLeptons.size() << std::endl;
    if (GenJetsWithoutLeptons.size() < 2)
    {
        return false;
    }

    float invMassLeadingJet = (GenJetsWithoutLeptons[0] + GenJetsWithoutLeptons[1]).M();
    if (invMassLeadingJet > minLeadingJetsInvMass)
        return true;
    else
        return false;
}
