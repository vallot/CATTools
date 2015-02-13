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
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "TH1F.h"
#include "TH2F.h"

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
  edm::EDGetTokenT<reco::GenParticleCollection> pseudotopToken_;

  double jetMinPt_, jetMaxEta_;
  double leptonMinPt_, leptonMaxEta_;
  TH2F* hL1PtParVsPseu_;
  TH2F* hL2PtParVsPseu_;
  TH2F* hB1PtParVsPseu_;
  TH2F* hB2PtParVsPseu_;
  
  TH1F* hPartonL1Pt_;
  TH1F* hPartonL1Eta_;
  TH1F* hPartonL1Phi_;
  TH1F* hPartonL2Pt_;
  TH1F* hPartonL2Eta_;
  TH1F* hPartonL2Phi_;

  TH1F* hPartonB1Pt_;
  TH1F* hPartonB1Eta_;
  TH1F* hPartonB1Phi_;

  TH1F* hPartonB2Pt_;
  TH1F* hPartonB2Eta_;
  TH1F* hPartonB2Phi_;

  TH1F* hPartonTop1Mass_;
  TH1F* hPartonTop2Mass_;

  TH1F* hPartonTTMass_;
  TH1F* hPartonTTDeltaPhi_;
  
  TH1F* hPartonLB1Mass_;
  TH1F* hPartonLB2Mass_;
  TH1F* hPartonDileptonMass_;

  TH1F* hPartonLBDeltaR_;
  
  TH1F* hPseudoL1Pt_;
  TH1F* hPseudoL1Eta_;
  TH1F* hPseudoL1Phi_;
  TH1F* hPseudoL2Pt_;
  TH1F* hPseudoL2Eta_;
  TH1F* hPseudoL2Phi_;

  TH1F* hPseudoB1Pt_;
  TH1F* hPseudoB1Eta_;
  TH1F* hPseudoB1Phi_;

  TH1F* hPseudoB2Pt_;
  TH1F* hPseudoB2Eta_;
  TH1F* hPseudoB2Phi_;

  TH1F* hPseudoTop1Mass_;
  TH1F* hPseudoTop2Mass_;

  TH1F* hPseudoTTMass_;
  TH1F* hPseudoTTDeltaPhi_;
  
  TH1F* hPseudoLB1Mass_;
  TH1F* hPseudoLB2Mass_;
  TH1F* hPseudoDileptonMass_;

  TH1F* hPseudoLBDeltaR_;
  
};

//using namespace cat;

CATGenLevelAnalysis::CATGenLevelAnalysis(const edm::ParameterSet& pset)
{
  partonsToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("partons"));
  pseudotopToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("pseudo"));
  channelToken_ = consumes<int>(pset.getParameter<edm::InputTag>("channel"));
  modesToken_ = consumes<std::vector<int> >(pset.getParameter<edm::InputTag>("modes"));

  leptonMinPt_ = pset.getParameter<double>("leptonMinPt");
  leptonMaxEta_ = pset.getParameter<double>("leptonMaxEta");


  jetMinPt_ = pset.getParameter<double>("jetMinPt");
  jetMaxEta_ = pset.getParameter<double>("jetMaxEta");

  edm::Service<TFileService> fs;
 
  hL1PtParVsPseu_ = fs->make<TH2F>("hL1PtParVsPseu", "lepton1 P_{T};parton p_{T} (GeV/c);pseudo p_{T} (GeV/c)", 100, 0, 150, 100, 0, 150);
  hL2PtParVsPseu_ = fs->make<TH2F>("hL2PtParVsPseu", "lepton2 P_{T};parton p_{T} (GeV/c);pseudo p_{T} (GeV/c)", 100, 0, 150, 100, 0, 150);
  hB1PtParVsPseu_ = fs->make<TH2F>("hB1PtParVsPseu", "b jet1 P_{T};parton p_{T} (GeV/c);pseudo p_{T} (GeV/c)", 100, 0, 200, 100, 0, 200);
  hB2PtParVsPseu_ = fs->make<TH2F>("hB2PtParVsPseu", "b jet2 P_{T};parton p_{T} (GeV/c);pseudo p_{T} (GeV/c)", 100, 0, 200, 100, 0, 200);
  
  TFileDirectory dirParton = fs->mkdir("parton", "parton");
  hPartonL1Pt_   = dirParton.make<TH1F>("hPartonL1Pt", "parton lepton 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPartonL1Eta_  = dirParton.make<TH1F>("hPartonL1Eta", "parton lepton 1 Eta;#eta;Events", 100, -3, 3);
  hPartonL1Phi_  = dirParton.make<TH1F>("hPartonL1Phi", "parton lepton 1 Phi;#phi;Events", 100, 0, acos(-1));

  hPartonL2Pt_   = dirParton.make<TH1F>("hPartonL2Pt", "parton lepton 2 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPartonL2Eta_  = dirParton.make<TH1F>("hPartonL2Eta", "parton lepton 2 Eta;#eta;Events", 100, -3, 3);
  hPartonL2Phi_  = dirParton.make<TH1F>("hPartonL2Phi", "parton lepton 2 Phi;#phi;Events", 100, 0, acos(-1));

  hPartonB1Pt_      = dirParton.make<TH1F>("hPartonB1Pt", "parton jet 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPartonB1Eta_     = dirParton.make<TH1F>("hPartonB1Eta", "parton jet 1 eta;#eta;Events", 100, -3, 3);
  hPartonB1Phi_     = dirParton.make<TH1F>("hPartonB1Phi", "parton jet 1 phi;#phi;Events", 100, 0, acos(-1));

  hPartonB2Pt_      = dirParton.make<TH1F>("hPartonB2Pt", "parton jet 2 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPartonB2Eta_     = dirParton.make<TH1F>("hPartonB2Eta", "parton jet 2 eta;#eta;Events", 100, -3, 3);
  hPartonB2Phi_     = dirParton.make<TH1F>("hPartonB2Phi", "parton jet 2 phi;#phi;Events", 100, 0, acos(-1));

  hPartonTop1Mass_      = dirParton.make<TH1F>("hPartonTop1Mass", "parton top 1 mass;Mass (GeV/c^{2};Events",100,100,300);
  hPartonTop2Mass_      = dirParton.make<TH1F>("hPartonTop2Mass", "parton top 2 mass;Mass (GeV/c^{2};Events",100,100,300);

  hPartonTTMass_      = dirParton.make<TH1F>("hPartonTTMass", "parton ttbar invariant mass;Mass (GeV/c^{2};Events",100,0,1000);
  hPartonTTDeltaPhi_  = dirParton.make<TH1F>("hPartonTTDeltaPhi", "parton ttbar #Delta#phi;#Delta#phi;Events",100,0,acos(-1));

  hPartonLB1Mass_  = dirParton.make<TH1F>("hPartonLB1Mass", "parton lepton 1 parton jet 1 invariant mass;Mass (GeV/c^{2});Events",100,0,200);
  hPartonLB2Mass_  = dirParton.make<TH1F>("hPartonLB2Mass", "parton lepton 2 parton jet 2 invariant mass;Mass (GeV/c^{2});Events",100,0,200);
  hPartonDileptonMass_  = dirParton.make<TH1F>("hPartonDileptonMass", "parton dilepton invariant mass;Mass (GeV/c^{2});Events",100,0,200);

  hPartonLBDeltaR_ = dirParton.make<TH1F>("hPartonLBDeltaR", "parton DeltaR;#DeltaR;Events", 100, 0, 5);

  TFileDirectory dirPseudo = fs->mkdir("pseudo", "pseudo");
  
  hPseudoL1Pt_   = dirPseudo.make<TH1F>("hPseudoL1Pt", "pseudo lepton 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPseudoL1Eta_  = dirPseudo.make<TH1F>("hPseudoL1Eta", "pseudo lepton 1 Eta;#eta;Events", 100, -3, 3);
  hPseudoL1Phi_  = dirPseudo.make<TH1F>("hPseudoL1Phi", "pseudo lepton 1 Phi;#phi;Events", 100, 0, acos(-1));

  hPseudoL2Pt_   = dirPseudo.make<TH1F>("hPseudoL2Pt", "pseudo lepton 2 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPseudoL2Eta_  = dirPseudo.make<TH1F>("hPseudoL2Eta", "pseudo lepton 2 Eta;#eta;Events", 100, -3, 3);
  hPseudoL2Phi_  = dirPseudo.make<TH1F>("hPseudoL2Phi", "pseudo lepton 2 Phi;#phi;Events", 100, 0, acos(-1));

  hPseudoB1Pt_      = dirPseudo.make<TH1F>("hPseudoB1Pt", "pseudo jet 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPseudoB1Eta_     = dirPseudo.make<TH1F>("hPseudoB1Eta", "pseudo jet 1 eta;#eta;Events", 100, -3, 3);
  hPseudoB1Phi_     = dirPseudo.make<TH1F>("hPseudoB1Phi", "pseudo jet 1 phi;#phi;Events", 100, 0, acos(-1));

  hPseudoB2Pt_      = dirPseudo.make<TH1F>("hPseudoB2Pt", "pseudo jet 2 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  hPseudoB2Eta_     = dirPseudo.make<TH1F>("hPseudoB2Eta", "pseudo jet 2 eta;#eta;Events", 100, -3, 3);
  hPseudoB2Phi_     = dirPseudo.make<TH1F>("hPseudoB2Phi", "pseudo jet 2 phi;#phi;Events", 100, 0, acos(-1));

  hPseudoTop1Mass_      = dirPseudo.make<TH1F>("hPseudoTop1Mass", "pseudo top 1 mass;Mass (GeV/c^{2};Events",100,100,300);
  hPseudoTop2Mass_      = dirPseudo.make<TH1F>("hPseudoTop2Mass", "pseudo top 2 mass;Mass (GeV/c^{2};Events",100,100,300);

  hPseudoTTMass_      = dirPseudo.make<TH1F>("hPseudoTTMass", "pseudo ttbar invariant mass;Mass (GeV/c^{2};Events",100,0,1000);
  hPseudoTTDeltaPhi_  = dirPseudo.make<TH1F>("hPseudoTTDeltaPhi", "pseudo ttbar #Delta#phi;#Delta#phi;Events",100,0,acos(-1));

  hPseudoLB1Mass_  = dirPseudo.make<TH1F>("hPseudoLB1Mass", "pseudo lepton 1 pseudo jet 1 invariant mass;Mass (GeV/c^{2});Events",100,0,200);
  hPseudoLB2Mass_  = dirPseudo.make<TH1F>("hPseudoLB2Mass", "pseudo lepton 2 pseudo jet 2 invariant mass;Mass (GeV/c^{2});Events",100,0,200);
  hPseudoDileptonMass_  = dirPseudo.make<TH1F>("hPseudoDileptonMass", "pseudo dilepton invariant mass;Mass (GeV/c^{2});Events",100,0,200);

  hPseudoLBDeltaR_ = dirPseudo.make<TH1F>("hPseudoLBDeltaR", "pseudo DeltaR;#DeltaR;Events", 100, 0, 5);
  
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

  edm::Handle<reco::GenParticleCollection> pseudoHandle;
  event.getByToken(pseudotopToken_, pseudoHandle);

  if (*channelHandle != 3) return;

  const reco::GenParticle* partonTop1 = 0, * partonTop2 = 0;
  for ( int i=0, n=partonsHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = partonsHandle->at(i);
    if ( p.pdgId() == 6 ) partonTop1 = &p;
    else if ( p.pdgId() == -6 ) partonTop2 = &p;
  }
  if ( !partonTop1 or !partonTop2 ) return;

  const reco::Candidate* partonW1 = partonTop1->daughter(0);
  const reco::Candidate* partonW2 = partonTop2->daughter(0);
  const reco::Candidate* partonB1 = partonTop1->daughter(1);
  const reco::Candidate* partonB2 = partonTop2->daughter(1);

  if (!partonW1 or !partonW2 or !partonB1 or !partonB2) return;

  const reco::Candidate* partonL1 = partonW1->daughter(0);
  const reco::Candidate* partonL2 = partonW2->daughter(0);
  if (!partonL1 or !partonL2) return;

  if ( abs(partonL1->pdgId()) == 15 and partonL1->numberOfDaughters() > 0 ) partonL1 = partonL1->daughter(0);
  if ( abs(partonL2->pdgId()) == 15 and partonL2->numberOfDaughters() > 0 ) partonL2 = partonL2->daughter(0);

  if ( partonL1->pt() < leptonMinPt_ or abs(partonL1->eta()) > leptonMaxEta_ ) return;
  if ( partonL2->pt() < leptonMinPt_ or abs(partonL2->eta()) > leptonMaxEta_ ) return;
  if ( partonB1->pt() < jetMinPt_ or abs(partonB1->eta()) > jetMaxEta_ ) return;
  if ( partonB2->pt() < jetMinPt_ or abs(partonB2->eta()) > jetMaxEta_ ) return;

  const reco::GenParticle* pseudoTop1 = 0, * pseudoTop2 = 0;
  for ( int i=0, n=pseudoHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = pseudoHandle->at(i);
    if ( p.pdgId() == 6 ) pseudoTop1 = &p;
    else if ( p.pdgId() == -6 ) pseudoTop2 = &p;
  }
  if ( !pseudoTop1 or !pseudoTop2 ) return;

  const reco::Candidate* pseudoW1 = pseudoTop1->daughter(0);
  const reco::Candidate* pseudoW2 = pseudoTop2->daughter(0);
  const reco::Candidate* pseudoB1 = pseudoTop1->daughter(1);
  const reco::Candidate* pseudoB2 = pseudoTop2->daughter(1);

  if (!pseudoW1 or !pseudoW2 or !pseudoB1 or !pseudoB2) return;

  const reco::Candidate* pseudoL1 = pseudoW1->daughter(0);
  const reco::Candidate* pseudoL2 = pseudoW2->daughter(0);
  if (!pseudoL1 or !pseudoL2) return;

  if ( abs(pseudoL1->pdgId()) == 15 and pseudoL1->numberOfDaughters() > 0 ) pseudoL1 = pseudoL1->daughter(0);
  if ( abs(pseudoL2->pdgId()) == 15 and pseudoL2->numberOfDaughters() > 0 ) pseudoL2 = pseudoL2->daughter(0);

  if ( pseudoL1->pt() < leptonMinPt_ or abs(pseudoL1->eta()) > leptonMaxEta_ ) return;
  if ( pseudoL2->pt() < leptonMinPt_ or abs(pseudoL2->eta()) > leptonMaxEta_ ) return;
  if ( pseudoB1->pt() < jetMinPt_ or abs(pseudoB1->eta()) > jetMaxEta_ ) return;
  if ( pseudoB2->pt() < jetMinPt_ or abs(pseudoB2->eta()) > jetMaxEta_ ) return;
  
  hL1PtParVsPseu_->Fill(partonL1->pt(),pseudoL1->pt());
  hL2PtParVsPseu_->Fill(partonL2->pt(),pseudoL2->pt());
  hB1PtParVsPseu_->Fill(partonB1->pt(),pseudoB1->pt());
  hB2PtParVsPseu_->Fill(partonB2->pt(),pseudoB2->pt());

  
  hPartonTop1Mass_->Fill(partonTop1->mass());
  hPartonTop2Mass_->Fill(partonTop2->mass());

  const double partonTTMass = (partonTop1->p4() + partonTop2->p4()).mass();
  const double partonTTDPhi = abs(deltaPhi(partonTop1->phi(), partonTop2->phi()));
  const double partonLBMass1 = (partonL1->p4() + partonB1->p4()).mass();
  const double partonLBMass2 = (partonL2->p4() + partonB2->p4()).mass();
  const double partonLLMass = (partonL1->p4() + partonL2->p4()).mass();

  hPartonL1Pt_->Fill(partonL1->pt());
  hPartonL1Eta_->Fill(partonL1->eta());
  hPartonL1Phi_->Fill(partonL1->phi());

  hPartonL2Pt_->Fill(partonL2->pt());
  hPartonL2Eta_->Fill(partonL2->eta());
  hPartonL2Phi_->Fill(partonL2->phi());

  hPartonB1Pt_->Fill(partonB1->pt());
  hPartonB1Eta_->Fill(partonB1->eta());
  hPartonB1Phi_->Fill(partonB1->phi());

  hPartonB2Pt_->Fill(partonB2->pt());
  hPartonB2Eta_->Fill(partonB2->eta());
  hPartonB2Phi_->Fill(partonB2->phi());

  hPartonTTMass_->Fill(partonTTMass);
  hPartonTTDeltaPhi_->Fill(partonTTDPhi);
  
  hPartonLB1Mass_->Fill(partonLBMass1);
  hPartonLB2Mass_->Fill(partonLBMass2);
  hPartonDileptonMass_->Fill(partonLLMass);

  hPartonLBDeltaR_->Fill(deltaR(partonL1->p4(),partonB1->p4()));    
  hPartonLBDeltaR_->Fill(deltaR(partonL2->p4(),partonB2->p4()));    
  
  hPseudoTop1Mass_->Fill(pseudoTop1->mass());
  hPseudoTop2Mass_->Fill(pseudoTop2->mass());

  const double pseudoTTMass = (pseudoTop1->p4() + pseudoTop2->p4()).mass();
  const double pseudoTTDPhi = abs(deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
  const double pseudoLBMass1 = (pseudoL1->p4() + pseudoB1->p4()).mass();
  const double pseudoLBMass2 = (pseudoL2->p4() + pseudoB2->p4()).mass();
  const double pseudoLLMass = (pseudoL1->p4() + pseudoL2->p4()).mass();

  hPseudoL1Pt_->Fill(pseudoL1->pt());
  hPseudoL1Eta_->Fill(pseudoL1->eta());
  hPseudoL1Phi_->Fill(pseudoL1->phi());

  hPseudoL2Pt_->Fill(pseudoL2->pt());
  hPseudoL2Eta_->Fill(pseudoL2->eta());
  hPseudoL2Phi_->Fill(pseudoL2->phi());

  hPseudoB1Pt_->Fill(pseudoB1->pt());
  hPseudoB1Eta_->Fill(pseudoB1->eta());
  hPseudoB1Phi_->Fill(pseudoB1->phi());

  hPseudoB2Pt_->Fill(pseudoB2->pt());
  hPseudoB2Eta_->Fill(pseudoB2->eta());
  hPseudoB2Phi_->Fill(pseudoB2->phi());

  hPseudoTTMass_->Fill(pseudoTTMass);
  hPseudoTTDeltaPhi_->Fill(pseudoTTDPhi);
  
  hPseudoLB1Mass_->Fill(pseudoLBMass1);
  hPseudoLB2Mass_->Fill(pseudoLBMass2);
  hPseudoDileptonMass_->Fill(pseudoLLMass);

  hPseudoLBDeltaR_->Fill(deltaR(pseudoL1->p4(),pseudoB1->p4()));    
  hPseudoLBDeltaR_->Fill(deltaR(pseudoL2->p4(),pseudoB2->p4()));    
  
}

DEFINE_FWK_MODULE(CATGenLevelAnalysis);

