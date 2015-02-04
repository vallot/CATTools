#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"

#include <iostream>

using namespace std;

class CATGenLevelAnalysis : public edm::EDAnalyzer
{
public:
  CATGenLevelAnalysis(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<int> channelToken_;
  edm::EDGetTokenT<std::vector<int> > modesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonsToken_;

  TH1F* hLepton1Pt_;
  TH1F* hLepton1Eta_;
  TH1F* hLepton1Phi_;
  TH1F* hLepton2Pt_;
  TH1F* hLepton2Eta_;
  TH1F* hLepton2Phi_;

  TH1F* hJet1Pt_;
  TH1F* hJet1Eta_;
  TH1F* hJet1Phi_;

  TH1F* hJet2Pt_;
  TH1F* hJet2Eta_;
  TH1F* hJet2Phi_;

  TH1F* hTTMass_;
  TH1F* hTTDeltaPhi_;
  
  TH1F* hLeptonJet1Mass_;
  TH1F* hLeptonJet2Mass_;
  TH1F* hDileptonMass_;

  TH1F* hLeptonJetDeltaR_;

};

//using namespace cat;

CATGenLevelAnalysis::CATGenLevelAnalysis(const edm::ParameterSet& pset)
{
  partonsToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("partons"));
  channelToken_ = consumes<int>(pset.getParameter<edm::InputTag>("channel"));
  modesToken_ = consumes<std::vector<int> >(pset.getParameter<edm::InputTag>("modes"));

  edm::Service<TFileService> fs;
  hLepton1Pt_   = fs->make<TH1F>("hLepton1Pt", "lepton 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hLepton1Eta_  = fs->make<TH1F>("hLepton1Eta", "lepton 1 Eta;#eta;Events", 100, -3, 3);
  hLepton1Phi_  = fs->make<TH1F>("hLepton1Phi", "lepton 1 Phi;#phi;Events", 100, 0, acos(-1));

  hLepton2Pt_   = fs->make<TH1F>("hLepton2Pt", "lepton 2 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hLepton2Eta_  = fs->make<TH1F>("hLepton2Eta", "lepton 2 Eta;#eta;Events", 100, -3, 3);
  hLepton2Phi_  = fs->make<TH1F>("hLepton2Phi", "lepton 2 Phi;#phi;Events", 100, 0, acos(-1));

  hJet1Pt_      = fs->make<TH1F>("hJet1Pt", "jet 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hJet1Eta_     = fs->make<TH1F>("hJet1Eta", "jet 1 eta;#eta;Events", 100, -3, 3);
  hJet1Phi_     = fs->make<TH1F>("hJet1Phi", "jet 1 phi;#phi;Events", 100, 0, acos(-1));

  hJet2Pt_      = fs->make<TH1F>("hJet2Pt", "jet 2 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hJet2Eta_     = fs->make<TH1F>("hJet2Eta", "jet 2 eta;#eta;Events", 100, -3, 3);
  hJet2Phi_     = fs->make<TH1F>("hJet2Phi", "jet 2 phi;#phi;Events", 100, 0, acos(-1));

  hTTMass_      = fs->make<TH1F>("hTTMass", "ttbar invariant mass;Mass (GeV/c^{2};Events",100,0,1000);
  hTTDeltaPhi_  = fs->make<TH1F>("hTTDeltaPhi", "ttbar #Delta#phi;#Delta#phi;Events",100,0,acos(-1));

  hLeptonJet1Mass_  = fs->make<TH1F>("hLeptonJet1Mass", "lepton 1 jet 1 invariant mass;Mass (GeV/c^{2});Events",100,0,200);
  hLeptonJet2Mass_  = fs->make<TH1F>("hLeptonJet2Mass", "lepton 2 jet 2 invariant mass;Mass (GeV/c^{2});Events",100,0,200);
  hDileptonMass_  = fs->make<TH1F>("hDileptonMass", "dilepton invariant mass;Mass (GeV/c^{2});Events",100,0,200);

  hLeptonJetDeltaR_ = fs->make<TH1F>("hLeptonJetDeltaR", "DeltaR;#DeltaR;Events", 100, 0, 5);
  // Do the same for the other histograms
}

void CATGenLevelAnalysis::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<int> channelHandle;
  event.getByToken(channelToken_, channelHandle);

  edm::Handle<std::vector<int> > modesHandle;
  event.getByToken(modesToken_, modesHandle);

  edm::Handle<reco::GenParticleCollection> partonsHandle;
  event.getByToken(partonsToken_, partonsHandle);

  if (*channelHandle != 3) return;
  cout << "CHANNEL = " << *channelHandle << endl;
  for ( int i=0, n=modesHandle->size(); i<n; ++i )
  {
    cout << "MODE" << i << " = " << modesHandle->at(i) << endl;
  }

  const reco::GenParticle* top1 = 0, * top2 = 0;
  for ( int i=0, n=partonsHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = partonsHandle->at(i);
    if ( p.pdgId() == 6 ) top1 = &p;
    else if ( p.pdgId() == -6 ) top2 = &p;
  }

  const reco::Candidate* w1 = top1->daughter(0);
  const reco::Candidate* w2 = top2->daughter(0);
  const reco::Candidate* b1 = top1->daughter(1);
  const reco::Candidate* b2 = top2->daughter(1);

  const reco::Candidate* l1 = w1->daughter(0);
  const reco::Candidate* l2 = w2->daughter(0);
  if ( abs(l1->pdgId()) == 15 and l1->numberOfDaughters() > 0 ) l1 = l1->daughter(0);
  if ( abs(l2->pdgId()) == 15 and l2->numberOfDaughters() > 0 ) l2 = l2->daughter(0);

  const double ttMass = (top1->p4() + top2->p4()).mass();
  const double ttDPhi = abs(deltaPhi(top1->phi(), top2->phi()));
  const double lbMass1 = (l1->p4() + b1->p4()).mass();
  const double lbMass2 = (l2->p4() + b2->p4()).mass();
  const double llMass = (l1->p4() + l2->p4()).mass();

  hLepton1Pt_->Fill(l1->pt());
  hLepton1Eta_->Fill(l1->eta());
  hLepton1Phi_->Fill(l1->phi());

  hLepton2Pt_->Fill(l2->pt());
  hLepton2Eta_->Fill(l2->eta());
  hLepton2Phi_->Fill(l2->phi());

  hJet1Pt_->Fill(b1->pt());
  hJet1Eta_->Fill(b1->eta());
  hJet1Phi_->Fill(b1->phi());

  hJet2Pt_->Fill(b2->pt());
  hJet2Eta_->Fill(b2->eta());
  hJet2Phi_->Fill(b2->phi());

  hTTMass_->Fill(ttMass);
  hTTDeltaPhi_->Fill(ttDPhi);
  
  hLeptonJet1Mass_->Fill(lbMass1);
  hLeptonJet2Mass_->Fill(lbMass2);
  hDileptonMass_->Fill(llMass);

  hLeptonJetDeltaR_->Fill(deltaR(l1->p4(),b1->p4()));    
  hLeptonJetDeltaR_->Fill(deltaR(l2->p4(),b2->p4()));    
  
}

DEFINE_FWK_MODULE(CATGenLevelAnalysis);

