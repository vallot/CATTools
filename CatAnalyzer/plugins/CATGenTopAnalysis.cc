#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "Math/GenVector/Boost.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;

class CATGenTopAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  CATGenTopAnalysis(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<int> channelToken_;
  edm::EDGetTokenT<std::vector<int> > modesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTopToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> pseudoTopToken_;

  const bool filterTaus_;

  typedef const reco::Candidate* CCandPtr;
  bool isAcceptedFullLept(CCandPtr l1, CCandPtr l2, CCandPtr b1, CCandPtr b2) const {
    if ( l1->pt() < 20 or std::abs(l1->eta()) > 2.4 ) return false;
    if ( l2->pt() < 20 or std::abs(l2->eta()) > 2.4 ) return false;
    if ( b1->pt() < 30 or std::abs(b1->eta()) > 2.4 ) return false;
    if ( b2->pt() < 30 or std::abs(b2->eta()) > 2.4 ) return false;
    return true;
  }
  bool isAcceptedSemiLept(CCandPtr lep, CCandPtr q1, CCandPtr q2, CCandPtr b1, CCandPtr b2) const {
    if ( lep->pt() < 33 or std::abs(lep->eta()) > 2.1 ) return false;
    if ( q1->pt() < 30 or std::abs(q1->eta()) > 2.4 ) return false;
    if ( q2->pt() < 30 or std::abs(q2->eta()) > 2.4 ) return false;
    if ( b1->pt() < 30 or std::abs(b1->eta()) > 2.4 ) return false;
    if ( b2->pt() < 30 or std::abs(b2->eta()) > 2.4 ) return false;
    return true;
  }

  typedef TH1F* H1;
  typedef TH2F* H2;

  enum {
    SL_topPt        , SL_topPtTtbarSys, SL_topY         ,
    SL_ttbarDelPhi  , SL_topPtLead    , SL_topPtSubLead ,
    SL_ttbarPt      , SL_ttbarY       , SL_ttbarMass    ,

    DL_topPt        , DL_topPtTtbarSys, DL_topY         ,
    DL_ttbarDelPhi  , DL_topPtLead    , DL_topPtSubLead ,
    DL_ttbarPt      , DL_ttbarY       , DL_ttbarMass    ,

    END
  };

  // 1D histograms
  std::vector<H1> hFulParton_; // Full phase space parton level
  std::vector<H1> hFidParton_; // Fiducial phase space parton level
  std::vector<H1> hComParton_; // Parton particle common phase space parton level
  std::vector<H1> hComPseudo_; // Parton particle common phase space particle level
  std::vector<H1> hPseudo_   ; // Particle level (fiducial cut by construction)
  std::vector<H1> hChPseudo_ ; // Particle level with decay parton level channel filter

  // Response matrices
  std::vector<H2> h2_, h2Com_;

  // Other plots for debugging
  // Decay channels
  H1 hFulParton_Channel_, hFidParton_Channel_;
  H1 hComParton_Channel_, hComPseudo_Channel_;
  H1 hPseudo_Channel_, hChPseudo_Channel_;

  // channel vs channel without looking at pseudotop phase space cut
  H2 h2DebugChannel_;
  H2 h2PseudoChannel_;
  H2 h2ComChannel_;
  H2 h2ChChannel_;

};

CATGenTopAnalysis::CATGenTopAnalysis(const edm::ParameterSet& pset):
  filterTaus_(pset.getParameter<bool>("filterTaus"))
{
  partonTopToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("partonTop"));
  pseudoTopToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("pseudoTop"));
  channelToken_ = consumes<int>(pset.getParameter<edm::InputTag>("channel"));
  modesToken_ = consumes<std::vector<int> >(pset.getParameter<edm::InputTag>("modes"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  const std::vector<std::vector<float> > bins = {
    {0, 60, 100, 150, 200, 260, 320, 400, 500}, // d15
    {0, 60, 100, 150, 200, 260, 320, 400, 500}, // d16
    {-2.5,-1.6,-1.2,-0.8,-0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.5}, // d17
    {0,2,2.75,3,3.15}, // d18
    {0,60,100,150,200,260,320,400,500}, //d19
    {0,60,100,150,200,260,320,400,500}, //d20
    {0.0,20,45,75,120,190,255}, //d21
    {-2.5,-1.3,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.3,2.5}, //d22
    {345,400,470,550,650,800,1100,1600}, //d23
    {0,65,125,200,290,400}, //d24
    {0,60,115,190,275,380,500}, //d25
    {-2.5,-1.6,-1,-0.5,0,0.5,1,1.6,2.5}, //d26
    {0,1.89,2.77,3.05,3.15}, //d27
    {0,75,130,200,290,400}, //d28
    {0,55,120,200,290,400}, //d29
    {0,30,80,170,300},
    {-2.5,-1.5,-1,-0.5,0,0.5,1,1.5,2.5},
    {340,380,470,620,820,1100,1600},
  };

  const std::vector<std::string> names = {
    "SL_topPt", "SL_topPtTtbarSys", "SL_topY", "SL_ttbarDelPhi",
    "SL_topPtLead", "SL_topPtSubLead",
    "SL_ttbarPt", "SL_ttbarY", "SL_ttbarMass",

    "DL_topPt", "DL_topPtTtbarSys", "DL_topY", "DL_ttbarDelPhi",
    "DL_topPtLead", "DL_topPtSubLead",
    "DL_ttbarPt", "DL_ttbarY", "DL_ttbarMass",
  };

  const std::vector<std::string> titles = {
    "top p_{T} (GeV)", "top p_{T} at CM frame (GeV)", "top rapidity", "#delta#phi(top1, top2)",
    "Leading top p_{T} (GeV)", "Subleading top p_{T} (GeV)",
    "t#bar{T} p_{T} (GeV)", "t#bar{t} rapidity (GeV)", "t#bar{t} mass (GeV)",

    "top p_{T} (GeV)", "top p_{T} at CM frame (GeV)", "top rapidity", "#delta#phi(top1, top2)",
    "Leading top p_{T} (GeV)", "Subleading top p_{T} (GeV)",
    "t#bar{T} p_{T} (GeV)", "t#bar{t} rapidity (GeV)", "t#bar{t} mass (GeV)",
  };

  assert(bins.size() == END);
  assert(names.size() == END);
  assert(titles.size() == END);

  auto dirFulParton = fs->mkdir("FullParton");
  auto dirFidParton = fs->mkdir("FiducialParton");
  auto dirComParton = fs->mkdir("CommonParton");
  auto dirComPseudo = fs->mkdir("CommonParticle");
  auto dirPseudo = fs->mkdir("Particle");
  auto dirChPseudo = fs->mkdir("ChFilteredParticle");

  for ( size_t i=0; i<END; ++i ) {
    const int nbins = bins[i].size()-1;
    const float* binPtr = &bins[i][0];

    const string name = names[i];
    const string title = string(name) + ";" + titles[i];
    hFulParton_.push_back(dirFulParton.make<TH1F>(name.c_str(), title.c_str(), nbins, binPtr));
    hFidParton_.push_back(dirFidParton.make<TH1F>(name.c_str(), title.c_str(), nbins, binPtr));
    hComParton_.push_back(dirComParton.make<TH1F>(name.c_str(), title.c_str(), nbins, binPtr));
    hComPseudo_.push_back(dirComPseudo.make<TH1F>(name.c_str(), title.c_str(), nbins, binPtr));
    hPseudo_.push_back(dirPseudo.make<TH1F>(name.c_str(), title.c_str(), nbins, binPtr));
    hChPseudo_.push_back(dirChPseudo.make<TH1F>(name.c_str(), title.c_str(), nbins, binPtr));

    const string name2 = "resp_"+names[i];
    const string title2 = names[i]+";Particle level "+titles[i]+";Parton level "+titles[i];
    h2_.push_back(dirPseudo.make<TH2F>(name2.c_str(), title2.c_str(), nbins, binPtr, nbins, binPtr));
    h2Com_.push_back(dirComPseudo.make<TH2F>(name2.c_str(), title2.c_str(), nbins, binPtr, nbins, binPtr));
  }

  assert(hFulParton_.size() == END);
  assert(hFidParton_.size() == END);
  assert(hComParton_.size() == END);
  assert(hComPseudo_.size() == END);
  assert(hPseudo_.size() == END);
  assert(hChPseudo_.size() == END);
  assert(h2_.size() == END);
  assert(h2Com_.size() == END);

  // Continue to debugging plots
  const std::vector<const char*> channelNames = {
    "Hadronic",
    "e+jet", "#mu+jet",
    "ee", "#mu#mu", "e#mu",
    "#tau hadronic+jet", "#tau#rightarrow l+jet", "#tau dilepton"
  };

  hFulParton_Channel_ = dirFulParton.make<TH1F>("channel", "channel", 9, 1, 10);
  hFidParton_Channel_ = dirFidParton.make<TH1F>("channel", "channel", 9, 1, 10);
  hComParton_Channel_ = dirComParton.make<TH1F>("channel", "channel", 9, 1, 10);
  hComPseudo_Channel_ = dirComPseudo.make<TH1F>("channel", "channel", 9, 1, 10);
  hPseudo_Channel_    = dirPseudo.make<TH1F>("channel", "channel", 9, 1, 10);
  hChPseudo_Channel_  = dirChPseudo.make<TH1F>("channel", "channel", 9, 1, 10);

  h2DebugChannel_  = fs->make<TH2F>("h2DebugChannel", "No filter debug;Parton;Pseudotop", 9, 1, 10, 9, 1, 10);
  h2PseudoChannel_ = fs->make<TH2F>("h2PseudoChannel", "Pseudotop phase space;Parton;Pseudotop", 9, 1, 10, 9, 1, 10);
  h2ComChannel_    = fs->make<TH2F>("h2ComChannel", "Common phase space;Parton;Pseudotop", 9, 1, 10, 9, 1, 10);
  h2ChChannel_     = fs->make<TH2F>("h2ChChannel", "Common phase space same channel;Parton;Pseudotop", 9, 1, 10, 9, 1, 10);

  std::vector<H2> hh = {h2DebugChannel_, h2PseudoChannel_, h2ComChannel_, h2ChChannel_};
  for ( int i=0, n=channelNames.size(); i<n; ++i ) {
    hFulParton_Channel_->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    hFidParton_Channel_->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    hComParton_Channel_->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    hComPseudo_Channel_->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    hPseudo_Channel_   ->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    hChPseudo_Channel_ ->GetXaxis()->SetBinLabel(i+1, channelNames[i]);

    h2DebugChannel_ ->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    h2PseudoChannel_->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    h2ComChannel_   ->GetXaxis()->SetBinLabel(i+1, channelNames[i]);
    h2ChChannel_    ->GetXaxis()->SetBinLabel(i+1, channelNames[i]);

    h2DebugChannel_ ->GetYaxis()->SetBinLabel(i+1, channelNames[i]);
    h2PseudoChannel_->GetYaxis()->SetBinLabel(i+1, channelNames[i]);
    h2ComChannel_   ->GetYaxis()->SetBinLabel(i+1, channelNames[i]);
    h2ChChannel_    ->GetYaxis()->SetBinLabel(i+1, channelNames[i]);
  }

}

void CATGenTopAnalysis::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<int> channelHandle;
  event.getByToken(channelToken_, channelHandle);
  const int channel = *channelHandle;
  if ( channel == 0 ) return;

  edm::Handle<std::vector<int> > modesHandle;
  event.getByToken(modesToken_, modesHandle);
  if ( modesHandle->size() != 2 ) return; // this should not happen if parton top module was not crashed
  const int mode1 = modesHandle->at(0), mode2 = modesHandle->at(1);
  if ( filterTaus_ and (mode1 >= 4 or mode2 >= 4) ) return;

  edm::Handle<reco::GenParticleCollection> partonTopHandle;
  event.getByToken(partonTopToken_, partonTopHandle);
  if ( partonTopHandle->empty() ) return;

  // Find all parton level top decay chain
  // Get Top quark pairs first
  const auto partonTop1 = &partonTopHandle->at(0);
  const auto partonTop2 = &partonTopHandle->at(1);
  if ( !partonTop1 or !partonTop2 ) return;
  // Get W and b quarks
  // Prdering is fixed from the PartonTopProducer, W first
  const auto partonW1 = partonTop1->daughter(0);
  const auto partonB1 = partonTop1->daughter(1);
  const auto partonW2 = partonTop2->daughter(0);
  const auto partonB2 = partonTop2->daughter(1);
  if ( !partonW1 or !partonB2 or !partonB1 or !partonB2 ) return;
  // Get W daughters
  // Ordering is fixed from the PartonTopProducer, lepton first.
  auto partonW11 = partonW1->daughter(0); // non-const to replace with its daughter if exists
  const auto partonW12 = partonW1->daughter(1);
  auto partonW21 = partonW2->daughter(0); // non-const to replace with its daughter if exists
  const auto partonW22 = partonW2->daughter(1);
  if ( !partonW11 or !partonW12 or !partonW21 or !partonW22 ) return;
  // Continue to daughter for leptonically decaying taus
  if ( abs(partonW11->pdgId()) == 15 and
       partonW11->numberOfDaughters() > 0 and
       abs(partonW11->daughter(0)->pdgId()) > 10 ) partonW11 = partonW11->daughter(0);
  if ( abs(partonW21->pdgId()) == 15 and
       partonW21->numberOfDaughters() > 0 and
       abs(partonW21->daughter(0)->pdgId()) > 10 ) partonW21 = partonW21->daughter(0);
  const auto partonTT = partonTop1->p4()+partonTop2->p4();
  const double partonTopPtAtCM = ROOT::Math::Boost(partonTT.BoostToCM())(partonTop1->p4()).pt();

  // Determine channel informations
  int partonTopCh = 0;
  if ( channel == 1 ) { // Full hadronic including taus
    if ( mode1 == 4 or mode2 == 4 ) partonTopCh = 7; // hadronic with tau
    else partonTopCh = 1;
  }
  else if ( channel == 2 ) { // semilepton channels
    if ( mode1 == 3 or mode2 == 3 ) partonTopCh = 2; // e+jets no tau
    else if ( mode1 == 2 or mode2 == 2 ) partonTopCh = 3; // mu+jets no tau
    else partonTopCh = 8; // any leptons from tau decay
  }
  else { // full leptonic channels
    if ( mode1 >= 4 or mode2 >= 4 ) partonTopCh = 9; // dilepton anything includes tau decay
    else if ( mode1 == 3 and mode2 == 3 ) partonTopCh = 4; // ee channel no tau
    else if ( mode1 == 2 and mode2 == 2 ) partonTopCh = 5; // mumu cnannel no tau
    else partonTopCh = 6; // emu channel no tau
  }
  hFulParton_Channel_->Fill(partonTopCh);

  // Fill parton top plots
  if ( partonTopCh == 2 or partonTopCh == 3 ) {
    hFulParton_[SL_topPt]->Fill(partonTop1->pt());
    hFulParton_[SL_topPt]->Fill(partonTop2->pt());
    hFulParton_[SL_topY]->Fill(partonTop1->p4().Rapidity());
    hFulParton_[SL_topY]->Fill(partonTop2->p4().Rapidity());
    hFulParton_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
    hFulParton_[SL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()));
    hFulParton_[SL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()));
    hFulParton_[SL_ttbarPt]->Fill(partonTT.pt());
    hFulParton_[SL_ttbarY]->Fill(partonTT.Rapidity());
    hFulParton_[SL_ttbarMass]->Fill(partonTT.mass());
    hFulParton_[SL_topPtTtbarSys]->Fill(partonTopPtAtCM);

    // Fill parton top plots in fiducial phase space
    if ( isAcceptedSemiLept(partonW11, partonW21, partonW22, partonB1, partonB2) ) {
      hFidParton_Channel_->Fill(partonTopCh);

      hFidParton_[SL_topPt]->Fill(partonTop1->pt());
      hFidParton_[SL_topPt]->Fill(partonTop2->pt());
      hFidParton_[SL_topY]->Fill(partonTop1->p4().Rapidity());
      hFidParton_[SL_topY]->Fill(partonTop2->p4().Rapidity());
      hFidParton_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
      hFidParton_[SL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()));
      hFidParton_[SL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()));
      hFidParton_[SL_ttbarPt]->Fill(partonTT.pt());
      hFidParton_[SL_ttbarY]->Fill(partonTT.Rapidity());
      hFidParton_[SL_ttbarMass]->Fill(partonTT.mass());
      hFidParton_[SL_topPtTtbarSys]->Fill(partonTopPtAtCM);
    }
  }
  else if ( partonTopCh == 4 or partonTopCh == 5 or partonTopCh == 6 ) {
    hFulParton_[DL_topPt]->Fill(partonTop1->pt());
    hFulParton_[DL_topPt]->Fill(partonTop2->pt());
    hFulParton_[DL_topY]->Fill(partonTop1->p4().Rapidity());
    hFulParton_[DL_topY]->Fill(partonTop2->p4().Rapidity());
    hFulParton_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
    hFulParton_[DL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()));
    hFulParton_[DL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()));
    hFulParton_[DL_ttbarPt]->Fill(partonTT.pt());
    hFulParton_[DL_ttbarY]->Fill(partonTT.Rapidity());
    hFulParton_[DL_ttbarMass]->Fill(partonTT.mass());
    hFulParton_[DL_topPtTtbarSys]->Fill(partonTopPtAtCM);

    if ( isAcceptedFullLept(partonW11, partonW21, partonB1, partonB2) ) {
      hFidParton_Channel_->Fill(partonTopCh);

      hFidParton_[DL_topPt]->Fill(partonTop1->pt());
      hFidParton_[DL_topPt]->Fill(partonTop2->pt());
      hFidParton_[DL_topY]->Fill(partonTop1->p4().Rapidity());
      hFidParton_[DL_topY]->Fill(partonTop2->p4().Rapidity());
      hFidParton_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
      hFidParton_[DL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()));
      hFidParton_[DL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()));
      hFidParton_[DL_ttbarPt]->Fill(partonTT.pt());
      hFidParton_[DL_ttbarY]->Fill(partonTT.Rapidity());
      hFidParton_[DL_ttbarMass]->Fill(partonTT.mass());
      hFidParton_[DL_topPtTtbarSys]->Fill(partonTopPtAtCM);
    }
  }

  edm::Handle<reco::GenParticleCollection> pseudoTopHandle;
  event.getByToken(pseudoTopToken_, pseudoTopHandle);
  if ( pseudoTopHandle->empty() )
  {
    // Fill parton top only (but no pseudoTop) plots
    return;
  }

  // Find all particle level top decay chain
  // Get Top quark pairs first
  const auto pseudoTop1 = &pseudoTopHandle->at(0);
  const auto pseudoTop2 = &pseudoTopHandle->at(1);

  // Get W and b quarks
  // Ordering is fixed, W first
  if ( !pseudoTop1 or !pseudoTop2 ) return;
  const auto pseudoW1 = pseudoTop1->daughter(0);
  const auto pseudoB1 = pseudoTop1->daughter(1);
  const auto pseudoW2 = pseudoTop2->daughter(0);
  const auto pseudoB2 = pseudoTop2->daughter(1);
  if ( !pseudoW1 or !pseudoW2 or !pseudoB1 or !pseudoB2 ) return;
  // Get W daughters
  // Ordering is fixed from the PartonTopProducer, lepton first.
  // There's no tau in PseudoTopProducer
  const auto pseudoW11 = pseudoW1->daughter(0);
  const auto pseudoW12 = pseudoW1->daughter(1);
  const auto pseudoW21 = pseudoW2->daughter(0);
  const auto pseudoW22 = pseudoW2->daughter(1);
  if ( !pseudoW11 or !pseudoW12 or !pseudoW21 or !pseudoW22 ) return;

  // Fill channel informations
  const int pseudoW1DauId = abs(pseudoW11->pdgId());
  const int pseudoW2DauId = abs(pseudoW21->pdgId());
  int pseudoTopCh = 0;
  if ( pseudoW1DauId < 10 and pseudoW2DauId < 10 ) pseudoTopCh = 1; // Full hadronic
  else if ( pseudoW1DauId > 10 and pseudoW2DauId > 10 ) { // dilepton
    switch ( pseudoW1DauId+pseudoW2DauId ) {
      case 22: pseudoTopCh = 4; break; // ee channel, 11+11
      case 26: pseudoTopCh = 5; break; // mumu channel, 13+13
      default: pseudoTopCh = 6; // others, emu channel
    }
  }
  else { // semilepton channel
    if ( pseudoW1DauId == 11 or pseudoW2DauId == 11 ) pseudoTopCh = 2; // e+jet
    else pseudoTopCh = 3; // mu+jet
  }
  hPseudo_Channel_->Fill(pseudoTopCh);
  h2DebugChannel_->Fill(partonTopCh, pseudoTopCh);

  const auto pseudoTT = pseudoTop1->p4()+pseudoTop2->p4();
  const double pseudoTopPtAtCM = ROOT::Math::Boost(pseudoTT.BoostToCM())(pseudoTop1->p4()).pt();

  // Fill pseudo top plots
  if ( pseudoTopCh == 1 ) { // Full hadronic in pseudoTop
    h2PseudoChannel_->Fill(partonTopCh, pseudoTopCh);
    //h2ComChannel_->Fill(partonTopCh, pseudoTopCh);
    if ( channel == 1 ) {
      hChPseudo_Channel_->Fill(pseudoTopCh);
      h2ChChannel_->Fill(partonTopCh, pseudoTopCh);
    }
  }
  else if ( (pseudoTopCh == 2 or pseudoTopCh == 3) and
            isAcceptedSemiLept(pseudoW11, pseudoW21, pseudoW22, pseudoB1, pseudoB2) ) {
    h2PseudoChannel_->Fill(partonTopCh, pseudoTopCh);

    // Additional acceptance cut for L+J channel
    hPseudo_[SL_topPt]->Fill(pseudoTop1->pt());
    hPseudo_[SL_topPt]->Fill(pseudoTop2->pt());
    hPseudo_[SL_topY]->Fill(pseudoTop1->p4().Rapidity());
    hPseudo_[SL_topY]->Fill(pseudoTop2->p4().Rapidity());
    hPseudo_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
    hPseudo_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()));
    hPseudo_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()));
    hPseudo_[SL_ttbarPt]->Fill(pseudoTT.pt());
    hPseudo_[SL_ttbarY]->Fill(pseudoTT.Rapidity());
    hPseudo_[SL_ttbarMass]->Fill(pseudoTT.mass());
    hPseudo_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM);

    // Fill response matrix no matter what parton level object acceptance is
    h2_[SL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt());
    h2_[SL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt());
    h2_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity());
    h2_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity());
    h2_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
    h2_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()));
    h2_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()));
    h2_[SL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt());
    h2_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity());
    h2_[SL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass());
    h2_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM);

    // Fill pseudo top plots within parton level acceptance cut
    if ( channel == 2 ) {
      hChPseudo_Channel_->Fill(pseudoTopCh);
      h2ChChannel_->Fill(partonTopCh, pseudoTopCh);

      hChPseudo_[SL_topPt]->Fill(pseudoTop1->pt());
      hChPseudo_[SL_topPt]->Fill(pseudoTop2->pt());
      hChPseudo_[SL_topY]->Fill(pseudoTop1->p4().Rapidity());
      hChPseudo_[SL_topY]->Fill(pseudoTop2->p4().Rapidity());
      hChPseudo_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
      hChPseudo_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()));
      hChPseudo_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()));
      hChPseudo_[SL_ttbarPt]->Fill(pseudoTT.pt());
      hChPseudo_[SL_ttbarY]->Fill(pseudoTT.Rapidity());
      hChPseudo_[SL_ttbarMass]->Fill(pseudoTT.mass());
      hChPseudo_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM);

      if ( isAcceptedSemiLept(partonW11, partonW21, partonW22, partonB1, partonB2) ) {
        hComParton_Channel_->Fill(partonTopCh);
        hComPseudo_Channel_->Fill(pseudoTopCh);
        h2ComChannel_->Fill(partonTopCh, pseudoTopCh);

        hComParton_[SL_topPt]->Fill(partonTop1->pt());
        hComParton_[SL_topPt]->Fill(partonTop2->pt());
        hComParton_[SL_topY]->Fill(partonTop1->p4().Rapidity());
        hComParton_[SL_topY]->Fill(partonTop2->p4().Rapidity());
        hComParton_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
        hComParton_[SL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()));
        hComParton_[SL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()));
        hComParton_[SL_ttbarPt]->Fill(partonTT.pt());
        hComParton_[SL_ttbarY]->Fill(partonTT.Rapidity());
        hComParton_[SL_ttbarMass]->Fill(partonTT.mass());
        hComParton_[SL_topPtTtbarSys]->Fill(partonTopPtAtCM);

        hComPseudo_[SL_topPt]->Fill(pseudoTop1->pt());
        hComPseudo_[SL_topPt]->Fill(pseudoTop2->pt());
        hComPseudo_[SL_topY]->Fill(pseudoTop1->p4().Rapidity());
        hComPseudo_[SL_topY]->Fill(pseudoTop2->p4().Rapidity());
        hComPseudo_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
        hComPseudo_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()));
        hComPseudo_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()));
        hComPseudo_[SL_ttbarPt]->Fill(pseudoTT.pt());
        hComPseudo_[SL_ttbarY]->Fill(pseudoTT.Rapidity());
        hComPseudo_[SL_ttbarMass]->Fill(pseudoTT.mass());
        hComPseudo_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM);

        // Fill response matrix no matter what parton level object acceptance is
        h2Com_[SL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt());
        h2Com_[SL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt());
        h2Com_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity());
        h2Com_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity());
        h2Com_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
        h2Com_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()));
        h2Com_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()));
        h2Com_[SL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt());
        h2Com_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity());
        h2Com_[SL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass());
        h2Com_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM);
      }
    }
  }
  else if ( (pseudoTopCh == 4 or pseudoTopCh == 5 or pseudoTopCh == 6 ) and
            isAcceptedFullLept(pseudoW11, pseudoW21, pseudoB1, pseudoB2) ) {
    h2PseudoChannel_->Fill(partonTopCh, pseudoTopCh);

    hPseudo_[DL_topPt]->Fill(pseudoTop1->pt());
    hPseudo_[DL_topPt]->Fill(pseudoTop2->pt());
    hPseudo_[DL_topY]->Fill(pseudoTop1->p4().Rapidity());
    hPseudo_[DL_topY]->Fill(pseudoTop2->p4().Rapidity());
    hPseudo_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
    hPseudo_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()));
    hPseudo_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()));
    hPseudo_[DL_ttbarPt]->Fill(pseudoTT.pt());
    hPseudo_[DL_ttbarY]->Fill(pseudoTT.Rapidity());
    hPseudo_[DL_ttbarMass]->Fill(pseudoTT.mass());
    hPseudo_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM);

    // Fill response matrix no matter what parton level object acceptance is
    h2_[DL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt());
    h2_[DL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt());
    h2_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity());
    h2_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity());
    h2_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
    h2_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()));
    h2_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()));
    h2_[DL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt());
    h2_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity());
    h2_[DL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass());
    h2_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM);

    if ( channel == 3 ) {
      hChPseudo_Channel_->Fill(pseudoTopCh);
      h2ChChannel_->Fill(partonTopCh, pseudoTopCh);

      hChPseudo_[DL_topPt]->Fill(pseudoTop1->pt());
      hChPseudo_[DL_topPt]->Fill(pseudoTop2->pt());
      hChPseudo_[DL_topY]->Fill(pseudoTop1->p4().Rapidity());
      hChPseudo_[DL_topY]->Fill(pseudoTop2->p4().Rapidity());
      hChPseudo_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
      hChPseudo_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()));
      hChPseudo_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()));
      hChPseudo_[DL_ttbarPt]->Fill(pseudoTT.pt());
      hChPseudo_[DL_ttbarY]->Fill(pseudoTT.Rapidity());
      hChPseudo_[DL_ttbarMass]->Fill(pseudoTT.mass());
      hChPseudo_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM);

      if ( isAcceptedFullLept(partonW11, partonW21, partonB1, partonB2) ) {
        hComParton_Channel_->Fill(partonTopCh);
        hComPseudo_Channel_->Fill(pseudoTopCh);
        h2ComChannel_->Fill(partonTopCh, pseudoTopCh);

        hComParton_[DL_topPt]->Fill(partonTop1->pt());
        hComParton_[DL_topPt]->Fill(partonTop2->pt());
        hComParton_[DL_topY]->Fill(partonTop1->p4().Rapidity());
        hComParton_[DL_topY]->Fill(partonTop2->p4().Rapidity());
        hComParton_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
        hComParton_[DL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()));
        hComParton_[DL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()));
        hComParton_[DL_ttbarPt]->Fill(partonTT.pt());
        hComParton_[DL_ttbarY]->Fill(partonTT.Rapidity());
        hComParton_[DL_ttbarMass]->Fill(partonTT.mass());
        hComParton_[DL_topPtTtbarSys]->Fill(partonTopPtAtCM);

        hComPseudo_[DL_topPt]->Fill(pseudoTop1->pt());
        hComPseudo_[DL_topPt]->Fill(pseudoTop2->pt());
        hComPseudo_[DL_topY]->Fill(pseudoTop1->p4().Rapidity());
        hComPseudo_[DL_topY]->Fill(pseudoTop2->p4().Rapidity());
        hComPseudo_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()));
        hComPseudo_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()));
        hComPseudo_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()));
        hComPseudo_[DL_ttbarPt]->Fill(pseudoTT.pt());
        hComPseudo_[DL_ttbarY]->Fill(pseudoTT.Rapidity());
        hComPseudo_[DL_ttbarMass]->Fill(pseudoTT.mass());
        hComPseudo_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM);

        // Fill response matrix no matter what parton level object acceptance is
        h2Com_[DL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt());
        h2Com_[DL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt());
        h2Com_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity());
        h2Com_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity());
        h2Com_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()));
        h2Com_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()));
        h2Com_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()));
        h2Com_[DL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt());
        h2Com_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity());
        h2Com_[DL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass());
        h2Com_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM);
      }
    }
  }

}

DEFINE_FWK_MODULE(CATGenTopAnalysis);

