#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;

class CATGenLeptonAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  CATGenLeptonAnalysis(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  TH1F* hZMass_, * hZDEta_, * hZDPhi_;
  TH1F* hZPt_, * hZEta_, * hZPhi_;
  TH1F* hL1Pt_, * hL1Eta_, * hL1Phi_;
  TH1F* hL2Pt_, * hL2Eta_, * hL2Phi_;
  TH1F* hLPPt_, * hLPEta_, * hLPPhi_;
  TH1F* hLMPt_, * hLMEta_, * hLMPhi_;

  TH1F* hNWZMass_, * hNWZDEta_, * hNWZDPhi_;
  TH1F* hNWZPt_, * hNWZEta_, * hNWZPhi_;
  TH1F* hNWL1Pt_, * hNWL1Eta_, * hNWL1Phi_;
  TH1F* hNWL2Pt_, * hNWL2Eta_, * hNWL2Phi_;
  TH1F* hNWLPPt_, * hNWLPEta_, * hNWLPPhi_;
  TH1F* hNWLMPt_, * hNWLMEta_, * hNWLMPhi_;

  TH2F* h2NuE_NuE_;
  TH2F* h2NWNuE_NuE_;

};

CATGenLeptonAnalysis::CATGenLeptonAnalysis(const edm::ParameterSet& pset)
{
  genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("src"));
  genWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("weight"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  auto dirDefault = fs->mkdir("default");
  hZMass_ = dirDefault.make<TH1F>("hZMass", "hZMass;M(l^{+}l^{-}) (GeV);Events per 1GeV", 1000, 0, 1000);
  hZDEta_ = dirDefault.make<TH1F>("hZDEta", "hZDEta;#delta#eta(l^{+},l^{-});Events per 0.05", 100, 0, 5);
  hZDPhi_ = dirDefault.make<TH1F>("hZDPhi", "hZDPhi;#delta#phi(l^{+},l^{-}) (radian);Events per 0.01", 315, 0, 3.15);

  hZPt_ = dirDefault.make<TH1F>("hZPt", "hZPt;p_{T}(l^{+}l^{-}) (GeV);Events per 1GeV", 1000, 0, 1000);
  hZEta_ = dirDefault.make<TH1F>("hZEta", "hZEta;#eta(l^{+}l^{-});Events per 0.05", 200, -5, 5);
  hZPhi_ = dirDefault.make<TH1F>("hZPhi", "hZPhi;#phi(l^{+}l^{-}) (radian);Events per 0.01", 2*315, -3.15, 3.15);

  hL1Pt_ = dirDefault.make<TH1F>("hL1Pt", "hL1Pt;p_{T}(l, 1st. lead) (GeV);Events per 1GeV", 1000, 0, 1000);
  hL2Pt_ = dirDefault.make<TH1F>("hL2Pt", "hL2Pt;p_{T}(l, 2nd. lead) (GeV);Events per 1GeV", 1000, 0, 1000);
  hLPPt_ = dirDefault.make<TH1F>("hLPPt", "hL+Pt;p_{T}(l^{+}) (GeV);Events per 1GeV", 1000, 0, 1000);
  hLMPt_ = dirDefault.make<TH1F>("hLMPt", "hL-Pt;p_{T}(l^{-}) (GeV);Events per 1GeV", 1000, 0, 1000);

  hL1Eta_ = dirDefault.make<TH1F>("hL1Eta", "hL1Eta;#eta(l, 1st. lead);Events per 0.05", 200, -5, 5);
  hL2Eta_ = dirDefault.make<TH1F>("hL2Eta", "hL2Eta;#eta(l, 2nd. lead);Events per 0.05", 200, -5, 5);
  hLPEta_ = dirDefault.make<TH1F>("hLPEta", "hL+Eta;#eta(l^{+});Events per 0.05", 200, -5, 5);
  hLMEta_ = dirDefault.make<TH1F>("hLMEta", "hL-Eta;#eta(l^{-});Events per 0.05", 200, -5, 5);

  hL1Phi_ = dirDefault.make<TH1F>("hL1Phi", "hL1Phi;#phi(l, 1st. lead) (radian);Events per 0.01", 2*315, -3.15, 3.15);
  hL2Phi_ = dirDefault.make<TH1F>("hL2Phi", "hL2Phi;#phi(l, 2nd. lead) (radian);Events per 0.01", 2*315, -3.15, 3.15);
  hLPPhi_ = dirDefault.make<TH1F>("hLPPhi", "hL+Phi;#phi(l^{+}) (radian);Events per 0.01", 2*315, -3.15, 3.15);
  hLMPhi_ = dirDefault.make<TH1F>("hLMPhi", "hL-Phi;#phi(l^{-}) (radian);Events per 0.01", 2*315, -3.15, 3.15);

  h2NuE_NuE_ = dirDefault.make<TH2F>("h2NuE_NuE", "h2NuE_NuE;Energy (#nu) (GeV);Energy (#bar{#nu}) (GeV)", 1000, 0, 1000, 1000, 0, 1000);

  auto dirNW = fs->mkdir("noweight");
  hNWZMass_ = dirNW.make<TH1F>("hZMass", "hZMass;M(l^{+}l^{-}) (GeV);Events per 1GeV", 1000, 0, 1000);
  hNWZDEta_ = dirNW.make<TH1F>("hZDEta", "hZDEta;#delta#eta(l^{+},l^{-});Events per 0.05", 100, 0, 5);
  hNWZDPhi_ = dirNW.make<TH1F>("hZDPhi", "hZDPhi;#delta#phi(l^{+},l^{-}) (radian);Events per 0.01", 315, 0, 3.15);

  hNWZPt_ = dirNW.make<TH1F>("hZPt", "hZPt;p_{T}(l^{+}l^{-}) (GeV);Events per 1GeV", 1000, 0, 1000);
  hNWZEta_ = dirNW.make<TH1F>("hZEta", "hZEta;#eta(l^{+}l^{-});Events per 0.05", 200, -5, 5);
  hNWZPhi_ = dirNW.make<TH1F>("hZPhi", "hZPhi;#phi(l^{+}l^{-}) (radian);Events per 0.01", 2*315, -3.15, 3.15);

  hNWL1Pt_ = dirNW.make<TH1F>("hL1Pt", "hL1Pt;p_{T}(l, 1st. lead) (GeV);Events per 1GeV", 1000, 0, 1000);
  hNWL2Pt_ = dirNW.make<TH1F>("hL2Pt", "hL2Pt;p_{T}(l, 2nd. lead) (GeV);Events per 1GeV", 1000, 0, 1000);
  hNWLPPt_ = dirNW.make<TH1F>("hLPPt", "hL+Pt;p_{T}(l^{+}) (GeV);Events per 1GeV", 1000, 0, 1000);
  hNWLMPt_ = dirNW.make<TH1F>("hLMPt", "hL-Pt;p_{T}(l^{-}) (GeV);Events per 1GeV", 1000, 0, 1000);

  hNWL1Eta_ = dirNW.make<TH1F>("hL1Eta", "hL1Eta;#eta(l, 1st. lead);Events per 0.05", 200, -5, 5);
  hNWL2Eta_ = dirNW.make<TH1F>("hL2Eta", "hL2Eta;#eta(l, 2nd. lead);Events per 0.05", 200, -5, 5);
  hNWLPEta_ = dirNW.make<TH1F>("hLPEta", "hL+Eta;#eta(l^{+});Events per 0.05", 200, -5, 5);
  hNWLMEta_ = dirNW.make<TH1F>("hLMEta", "hL-Eta;#eta(l^{-});Events per 0.05", 200, -5, 5);

  hNWL1Phi_ = dirNW.make<TH1F>("hL1Phi", "hL1Phi;#phi(l, 1st. lead) (radian);Events per 0.01", 2*315, -3.15, 3.15);
  hNWL2Phi_ = dirNW.make<TH1F>("hL2Phi", "hL2Phi;#phi(l, 2nd. lead) (radian);Events per 0.01", 2*315, -3.15, 3.15);
  hNWLPPhi_ = dirNW.make<TH1F>("hLPPhi", "hL+Phi;#phi(l^{+}) (radian);Events per 0.01", 2*315, -3.15, 3.15);
  hNWLMPhi_ = dirNW.make<TH1F>("hLMPhi", "hL-Phi;#phi(l^{-}) (radian);Events per 0.01", 2*315, -3.15, 3.15);

  h2NWNuE_NuE_ = dirDefault.make<TH2F>("h2NuE_NuE", "h2NuE_NuE;Energy (#nu) (GeV);Energy (#bar{#nu}) (GeV)", 1000, 0, 1000, 1000, 0, 1000);

}

void CATGenLeptonAnalysis::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<float> genWeightHandle;
  event.getByToken(genWeightToken_, genWeightHandle);
  const auto weight = *genWeightHandle;

  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByToken(genParticlesToken_, genParticlesHandle);

  // Collect generator level leptons, before radiation
  typedef const reco::GenParticle* GenParticlePtr;
  std::vector<GenParticlePtr> lepPs, lepMs;
  std::vector<GenParticlePtr> nus, nubars;
  for ( const auto& p : *genParticlesHandle )
  {
    if ( p.status() == 1 ) continue;
    const int id = p.pdgId();
    if ( p.numberOfMothers() > 0 and p.mother()->pdgId() == id ) continue;
    const int absId = abs(id);

    if ( absId == 11 or absId == 13 or absId == 15 )
    {
      if ( abs(p.mother()->pdgId()) == 15 ) continue;  //veto e/mu from tau.
      if ( p.charge() > 0 ) lepPs.push_back(&p);
      else lepMs.push_back(&p);
    }
    else if ( id ==  12 or id ==  14 ) nus.push_back(&p);
    else if ( id == -12 or id == -14 ) nubars.push_back(&p);
  }
  if ( lepPs.empty() or lepMs.empty() ) return;

  auto greaterByPtPtr = [](GenParticlePtr a, GenParticlePtr b){return a->pt() > b->pt();};
  std::nth_element(lepPs.begin(), lepPs.begin()+1, lepPs.end(), greaterByPtPtr);
  std::nth_element(lepMs.begin(), lepMs.begin()+1, lepMs.end(), greaterByPtPtr);

  const auto lepPLV = lepPs.at(0)->p4(), lepMLV = lepMs.at(0)->p4();
  const auto lep1LV = lepPLV.pt() > lepMLV.pt() ? lepPLV : lepMLV;
  const auto lep2LV = lepPLV.pt() > lepMLV.pt() ? lepMLV : lepPLV;
  const auto zLV = lepPLV + lepMLV;
  const double dEta = std::abs(lepPLV.eta()-lepMLV.eta());
  const double dPhi = std::abs(reco::deltaPhi(lepPLV.phi(), lepMLV.phi()));

  hZMass_->Fill(zLV.mass(), weight);
  hZDEta_->Fill(dEta, weight);
  hZDPhi_->Fill(dPhi, weight);

  hZPt_->Fill(zLV.pt(), weight);
  hZEta_->Fill(zLV.eta(), weight);
  hZPhi_->Fill(zLV.phi(), weight);

  hL1Pt_->Fill(lep1LV.pt(), weight);
  hL1Eta_->Fill(lep1LV.eta(), weight);
  hL1Phi_->Fill(lep1LV.phi(), weight);

  hL2Pt_->Fill(lep2LV.pt(), weight);
  hL2Eta_->Fill(lep2LV.eta(), weight);
  hL2Phi_->Fill(lep2LV.phi(), weight);

  hLPPt_->Fill(lepPLV.pt(), weight);
  hLPEta_->Fill(lepPLV.eta(), weight);
  hLPPhi_->Fill(lepPLV.phi(), weight);

  hLMPt_->Fill(lepMLV.pt(), weight);
  hLMEta_->Fill(lepMLV.eta(), weight);
  hLMPhi_->Fill(lepMLV.phi(), weight);

  hNWZMass_->Fill(zLV.mass());
  hNWZDEta_->Fill(dEta);
  hNWZDPhi_->Fill(dPhi);

  hNWZPt_->Fill(zLV.pt());
  hNWZEta_->Fill(zLV.eta());
  hNWZPhi_->Fill(zLV.phi());

  hNWL1Pt_->Fill(lep1LV.pt());
  hNWL1Eta_->Fill(lep1LV.eta());
  hNWL1Phi_->Fill(lep1LV.phi());

  hNWL2Pt_->Fill(lep2LV.pt());
  hNWL2Eta_->Fill(lep2LV.eta());
  hNWL2Phi_->Fill(lep2LV.phi());

  hNWLPPt_->Fill(lepPLV.pt());
  hNWLPEta_->Fill(lepPLV.eta());
  hNWLPPhi_->Fill(lepPLV.phi());

  hNWLMPt_->Fill(lepMLV.pt());
  hNWLMEta_->Fill(lepMLV.eta());
  hNWLMPhi_->Fill(lepMLV.phi());

  // Fill neutrinos
  if ( !nus.empty() and !nubars.empty() and
       abs(lepPs.at(0)->pdgId()) != 15 and abs(lepMs.at(0)->pdgId()) != 15 )
  {
    std::nth_element(nus.begin(), nus.begin()+1, nus.end(), greaterByPtPtr);
    std::nth_element(nubars.begin(), nubars.begin()+1, nubars.end(), greaterByPtPtr);

    h2NuE_NuE_->Fill(nus.at(0)->energy(), nubars.at(0)->energy(), weight);
    h2NWNuE_NuE_->Fill(nus.at(0)->energy(), nubars.at(0)->energy());
  }
}

DEFINE_FWK_MODULE(CATGenLeptonAnalysis);

