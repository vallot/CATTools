#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"

#include "Math/GenVector/Boost.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>

using namespace std;
using namespace cat;

class CATGenTopAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  CATGenTopAnalysis(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup&) override;

private:
  typedef std::vector<float> vfloat;

  edm::EDGetTokenT<float> weightToken_;
  edm::EDGetTokenT<vfloat> weightsToken_;
  int weightIndex_;

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

  typedef TH1D* H1;
  typedef TH2D* H2;

  enum H {
    SL_topPt        , SL_topPtTtbarSys, SL_topY         ,
    SL_ttbarDelPhi  , SL_topPtLead    , SL_topPtSubLead ,
    SL_ttbarPt      , SL_ttbarY       , SL_ttbarMass    ,

    DL_topPt        , DL_topPtTtbarSys, DL_topY         ,
    DL_ttbarDelPhi  , DL_topPtLead    , DL_topPtSubLead ,
    DL_ttbarPt      , DL_ttbarY       , DL_ttbarMass    ,

    END
  };

  // 1D histograms
  H1 hWeight_; // Weight distribution

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

  weightIndex_ = pset.getParameter<int>("weightIndex");
  if ( weightIndex_ < 0 ) weightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("weight"));
  else weightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("weight"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  typedef std::tuple<std::string, std::string, int, double, double> HDEF;
  std::map<H, HDEF> hdefs;
  hdefs[SL_topPt         ] = HDEF("SL_topPt", "top p_{T} (GeV)", 1000, 0, 1000);
  hdefs[SL_topPtTtbarSys ] = HDEF("SL_topPtTtbarSys", "top p_{T} at CM frame (GeV)", 1000, 0, 1000);
  hdefs[SL_topY          ] = HDEF("SL_topY", "top rapidity", 100, -5, 5);
  hdefs[SL_ttbarDelPhi   ] = HDEF("SL_ttbarDelPhi", "#delta#phi(top1, top2)", 100, 0, 3.15);
  hdefs[SL_topPtLead     ] = HDEF("SL_topPtLead", "Leading top p_{T} (GeV)", 1000, 0, 1000);
  hdefs[SL_topPtSubLead  ] = HDEF("SL_topPtSubLead", "Subleading top p_{T} (GeV)", 1000, 0, 1000);
  hdefs[SL_ttbarPt       ] = HDEF("SL_ttbarPt", "t#bar{t} p_{T} (GeV)", 1000, 0, 1000);
  hdefs[SL_ttbarY        ] = HDEF("SL_ttbarY", "t#bar{t} rapidity (GeV)", 100, -5, 5);
  hdefs[SL_ttbarMass     ] = HDEF("SL_ttbarMass", "t#bar{t} mass (GeV)", 5000, 0, 5000);

  hdefs[DL_topPt         ] = HDEF("DL_topPt", "top p_{T} (GeV)", 1000, 0, 1000);
  hdefs[DL_topPtTtbarSys ] = HDEF("DL_topPtTtbarSys", "top p_{T} at CM frame (GeV)", 1000, 0, 1000);
  hdefs[DL_topY          ] = HDEF("DL_topY", "top rapidity", 100, -5, 5);
  hdefs[DL_ttbarDelPhi   ] = HDEF("DL_ttbarDelPhi", "#delta#phi(top1, top2)", 100, 0, 3.15);
  hdefs[DL_topPtLead     ] = HDEF("DL_topPtLead", "Leading top p_{T} (GeV)", 1000, 0, 1000);
  hdefs[DL_topPtSubLead  ] = HDEF("DL_topPtSubLead", "Subleading top p_{T} (GeV)", 1000, 0, 1000);
  hdefs[DL_ttbarPt       ] = HDEF("DL_ttbarPt", "t#bar{t} p_{T} (GeV)", 1000, 0, 1000);
  hdefs[DL_ttbarY        ] = HDEF("DL_ttbarY", "t#bar{t} rapidity (GeV)", 100, -5, 5);
  hdefs[DL_ttbarMass     ] = HDEF("DL_ttbarMass", "t#bar{t} mass (GeV)", 5000, 0, 5000);

  assert(hdefs.size() == END);

  auto dirFulParton = fs->mkdir("FullParton");
  auto dirFidParton = fs->mkdir("FiducialParton");
  auto dirComParton = fs->mkdir("CommonParton");
  auto dirComPseudo = fs->mkdir("CommonParticle");
  auto dirPseudo = fs->mkdir("Particle");
  auto dirChPseudo = fs->mkdir("ChFilteredParticle");

  for ( auto x : hdefs ) {
    const auto name = std::get<0>(x.second);
    const auto title0 = std::get<1>(x.second);
    const auto title = name+";"+title0;
    const auto nbins = std::get<2>(x.second);
    const auto xlo = std::get<3>(x.second);
    const auto xhi = std::get<4>(x.second);

    hFulParton_.push_back(dirFulParton.make<TH1D>(name.c_str(), title.c_str(), nbins, xlo, xhi));
    hFidParton_.push_back(dirFidParton.make<TH1D>(name.c_str(), title.c_str(), nbins, xlo, xhi));
    hComParton_.push_back(dirComParton.make<TH1D>(name.c_str(), title.c_str(), nbins, xlo, xhi));
    hComPseudo_.push_back(dirComPseudo.make<TH1D>(name.c_str(), title.c_str(), nbins, xlo, xhi));
    hPseudo_.push_back(dirPseudo.make<TH1D>(name.c_str(), title.c_str(), nbins, xlo, xhi));
    hChPseudo_.push_back(dirChPseudo.make<TH1D>(name.c_str(), title.c_str(), nbins, xlo, xhi));

    const string name2 = "resp_"+name;
    const string title2 = name+";Particle level "+title0+";Parton level "+title0;
    h2_.push_back(dirPseudo.make<TH2D>(name2.c_str(), title2.c_str(), nbins, xlo, xhi, nbins, xlo, xhi));
    h2Com_.push_back(dirComPseudo.make<TH2D>(name2.c_str(), title2.c_str(), nbins, xlo, xhi, nbins, xlo, xhi));
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

  hFulParton_Channel_ = dirFulParton.make<TH1D>("channel", "channel", 9, 0, 9);
  hFidParton_Channel_ = dirFidParton.make<TH1D>("channel", "channel", 9, 0, 9);
  hComParton_Channel_ = dirComParton.make<TH1D>("channel", "channel", 9, 0, 9);
  hComPseudo_Channel_ = dirComPseudo.make<TH1D>("channel", "channel", 9, 0, 9);
  hPseudo_Channel_    = dirPseudo.make<TH1D>("channel", "channel", 9, 0, 9);
  hChPseudo_Channel_  = dirChPseudo.make<TH1D>("channel", "channel", 9, 0, 9);

  h2DebugChannel_  = fs->make<TH2D>("h2DebugChannel", "No filter debug;Parton;Pseudotop", 9, 0, 9, 9, 0, 9);
  h2PseudoChannel_ = fs->make<TH2D>("h2PseudoChannel", "Pseudotop phase space;Parton;Pseudotop", 9, 0, 9, 9, 0, 9);
  h2ComChannel_    = fs->make<TH2D>("h2ComChannel", "Common phase space;Parton;Pseudotop", 9, 0,9, 9, 0, 9);
  h2ChChannel_     = fs->make<TH2D>("h2ChChannel", "Common phase space same channel;Parton;Pseudotop", 9, 0, 9, 9, 0, 9);

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

  hWeight_ = fs->make<TH1D>("hWeight", "weights", 100, -5, 5);

}

void CATGenTopAnalysis::analyze(const edm::Event& event, const edm::EventSetup&)
{
  float weight = 1.;
  if ( weightIndex_ < 0 ) {
    edm::Handle<float> handle;
    event.getByToken(weightToken_, handle);
    weight = *handle;
  }
  else {
    edm::Handle<vfloat> handle;
    event.getByToken(weightsToken_, handle);
    weight = handle->at(weightIndex_);
  }

  hWeight_->Fill(weight);

  edm::Handle<int> channelHandle;
  event.getByToken(channelToken_, channelHandle);
  const int channel = *channelHandle;
  if ( channel == CH_NOTT ) return;

  edm::Handle<std::vector<int> > modesHandle;
  event.getByToken(modesToken_, modesHandle);
  if ( modesHandle->size() != 2 ) return; // this should not happen if parton top module was not crashed
  const int mode1 = modesHandle->at(0), mode2 = modesHandle->at(1);
  if ( filterTaus_ and (mode1 >= CH_TAU_HADRON or mode2 >= CH_TAU_HADRON) ) return;

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
  int partonTopCh = -1;
  if ( channel == CH_FULLHADRON ) { // Full hadronic including taus
    if ( mode1 == CH_TAU_HADRON or mode2 == CH_TAU_HADRON ) partonTopCh = 6; // hadronic with tau
    else partonTopCh = 0;
  }
  else if ( channel == CH_SEMILEPTON ) { // semilepton channels
    if ( mode1 == CH_ELECTRON or mode2 == CH_ELECTRON ) partonTopCh = 1; // e+jets no tau
    else if ( mode1 == CH_MUON or mode2 == CH_MUON ) partonTopCh = 2; // mu+jets no tau
    else partonTopCh = 7; // any leptons from tau decay
  }
  else { // full leptonic channels
    if ( mode1 > CH_TAU_HADRON or mode2 > CH_TAU_HADRON ) partonTopCh = 8; // dilepton anything includes tau decay
    else if ( mode1 == CH_ELECTRON and mode2 == CH_ELECTRON ) partonTopCh = 3; // ee channel no tau
    else if ( mode1 == CH_MUON and mode2 == CH_MUON ) partonTopCh = 4; // mumu cnannel no tau
    else partonTopCh = 5; // emu channel no tau
  }
  hFulParton_Channel_->Fill(partonTopCh, weight);

  // Fill parton top plots
  if ( partonTopCh == 1 or partonTopCh == 2 ) {
    hFulParton_[SL_topPt]->Fill(partonTop1->pt(), weight);
    hFulParton_[SL_topPt]->Fill(partonTop2->pt(), weight);
    hFulParton_[SL_topY]->Fill(partonTop1->p4().Rapidity(), weight);
    hFulParton_[SL_topY]->Fill(partonTop2->p4().Rapidity(), weight);
    hFulParton_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
    hFulParton_[SL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()), weight);
    hFulParton_[SL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()), weight);
    hFulParton_[SL_ttbarPt]->Fill(partonTT.pt(), weight);
    hFulParton_[SL_ttbarY]->Fill(partonTT.Rapidity(), weight);
    hFulParton_[SL_ttbarMass]->Fill(partonTT.mass(), weight);
    hFulParton_[SL_topPtTtbarSys]->Fill(partonTopPtAtCM, weight);

    // Fill parton top plots in fiducial phase space
    if ( isAcceptedSemiLept(partonW11, partonW21, partonW22, partonB1, partonB2) ) {
      hFidParton_Channel_->Fill(partonTopCh, weight);

      hFidParton_[SL_topPt]->Fill(partonTop1->pt(), weight);
      hFidParton_[SL_topPt]->Fill(partonTop2->pt(), weight);
      hFidParton_[SL_topY]->Fill(partonTop1->p4().Rapidity(), weight);
      hFidParton_[SL_topY]->Fill(partonTop2->p4().Rapidity(), weight);
      hFidParton_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
      hFidParton_[SL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()), weight);
      hFidParton_[SL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()), weight);
      hFidParton_[SL_ttbarPt]->Fill(partonTT.pt(), weight);
      hFidParton_[SL_ttbarY]->Fill(partonTT.Rapidity(), weight);
      hFidParton_[SL_ttbarMass]->Fill(partonTT.mass(), weight);
      hFidParton_[SL_topPtTtbarSys]->Fill(partonTopPtAtCM, weight);
    }
  }
  else if ( partonTopCh == 3 or partonTopCh == 4 or partonTopCh == 5 ) {
    hFulParton_[DL_topPt]->Fill(partonTop1->pt(), weight);
    hFulParton_[DL_topPt]->Fill(partonTop2->pt(), weight);
    hFulParton_[DL_topY]->Fill(partonTop1->p4().Rapidity(), weight);
    hFulParton_[DL_topY]->Fill(partonTop2->p4().Rapidity(), weight);
    hFulParton_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
    hFulParton_[DL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()), weight);
    hFulParton_[DL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()), weight);
    hFulParton_[DL_ttbarPt]->Fill(partonTT.pt(), weight);
    hFulParton_[DL_ttbarY]->Fill(partonTT.Rapidity(), weight);
    hFulParton_[DL_ttbarMass]->Fill(partonTT.mass(), weight);
    hFulParton_[DL_topPtTtbarSys]->Fill(partonTopPtAtCM, weight);

    if ( isAcceptedFullLept(partonW11, partonW21, partonB1, partonB2) ) {
      hFidParton_Channel_->Fill(partonTopCh, weight);

      hFidParton_[DL_topPt]->Fill(partonTop1->pt(), weight);
      hFidParton_[DL_topPt]->Fill(partonTop2->pt(), weight);
      hFidParton_[DL_topY]->Fill(partonTop1->p4().Rapidity(), weight);
      hFidParton_[DL_topY]->Fill(partonTop2->p4().Rapidity(), weight);
      hFidParton_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
      hFidParton_[DL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()), weight);
      hFidParton_[DL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()), weight);
      hFidParton_[DL_ttbarPt]->Fill(partonTT.pt(), weight);
      hFidParton_[DL_ttbarY]->Fill(partonTT.Rapidity(), weight);
      hFidParton_[DL_ttbarMass]->Fill(partonTT.mass(), weight);
      hFidParton_[DL_topPtTtbarSys]->Fill(partonTopPtAtCM, weight);
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
  int pseudoTopCh = -1;
  if ( pseudoW1DauId < 10 and pseudoW2DauId < 10 ) pseudoTopCh = 0; // Full hadronic
  else if ( pseudoW1DauId > 10 and pseudoW2DauId > 10 ) { // dilepton
    switch ( pseudoW1DauId+pseudoW2DauId ) {
      case 22: pseudoTopCh = 3; break; // ee channel, 11+11
      case 26: pseudoTopCh = 4; break; // mumu channel, 13+13
      default: pseudoTopCh = 5; // others, emu channel
    }
  }
  else { // semilepton channel
    if ( pseudoW1DauId == 11 or pseudoW2DauId == 11 ) pseudoTopCh = 1; // e+jet
    else pseudoTopCh = 2; // mu+jet
  }
  hPseudo_Channel_->Fill(pseudoTopCh, weight);
  h2DebugChannel_->Fill(partonTopCh, pseudoTopCh, weight);

  const auto pseudoTT = pseudoTop1->p4()+pseudoTop2->p4();
  const double pseudoTopPtAtCM = ROOT::Math::Boost(pseudoTT.BoostToCM())(pseudoTop1->p4()).pt();

  // Fill pseudo top plots
  if ( pseudoTopCh == 0 ) { // Full hadronic in pseudoTop
    h2PseudoChannel_->Fill(partonTopCh, pseudoTopCh, weight);
    //h2ComChannel_->Fill(partonTopCh, pseudoTopCh, weight);
    if ( channel == 0 ) {
      hChPseudo_Channel_->Fill(pseudoTopCh, weight);
      h2ChChannel_->Fill(partonTopCh, pseudoTopCh, weight);
    }
  }
  else if ( (pseudoTopCh == 1 or pseudoTopCh == 2) and
            isAcceptedSemiLept(pseudoW11, pseudoW21, pseudoW22, pseudoB1, pseudoB2) ) {
    h2PseudoChannel_->Fill(partonTopCh, pseudoTopCh, weight);

    // Additional acceptance cut for L+J channel
    hPseudo_[SL_topPt]->Fill(pseudoTop1->pt(), weight);
    hPseudo_[SL_topPt]->Fill(pseudoTop2->pt(), weight);
    hPseudo_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), weight);
    hPseudo_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), weight);
    hPseudo_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), weight);
    hPseudo_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), weight);
    hPseudo_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), weight);
    hPseudo_[SL_ttbarPt]->Fill(pseudoTT.pt(), weight);
    hPseudo_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), weight);
    hPseudo_[SL_ttbarMass]->Fill(pseudoTT.mass(), weight);
    hPseudo_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, weight);

    // Fill response matrix no matter what parton level object acceptance is
    h2_[SL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt(), weight);
    h2_[SL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt(), weight);
    h2_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
    h2_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
    h2_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
    h2_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()), weight);
    h2_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()), weight);
    h2_[SL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt(), weight);
    h2_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity(), weight);
    h2_[SL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass(), weight);
    h2_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM, weight);

    // Fill pseudo top plots within parton level acceptance cut
    if ( channel == CH_SEMILEPTON ) {
      hChPseudo_Channel_->Fill(pseudoTopCh, weight);
      h2ChChannel_->Fill(partonTopCh, pseudoTopCh, weight);

      hChPseudo_[SL_topPt]->Fill(pseudoTop1->pt(), weight);
      hChPseudo_[SL_topPt]->Fill(pseudoTop2->pt(), weight);
      hChPseudo_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), weight);
      hChPseudo_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), weight);
      hChPseudo_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), weight);
      hChPseudo_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), weight);
      hChPseudo_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), weight);
      hChPseudo_[SL_ttbarPt]->Fill(pseudoTT.pt(), weight);
      hChPseudo_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), weight);
      hChPseudo_[SL_ttbarMass]->Fill(pseudoTT.mass(), weight);
      hChPseudo_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, weight);

      if ( isAcceptedSemiLept(partonW11, partonW21, partonW22, partonB1, partonB2) ) {
        hComParton_Channel_->Fill(partonTopCh, weight);
        hComPseudo_Channel_->Fill(pseudoTopCh, weight);
        h2ComChannel_->Fill(partonTopCh, pseudoTopCh, weight);

        hComParton_[SL_topPt]->Fill(partonTop1->pt(), weight);
        hComParton_[SL_topPt]->Fill(partonTop2->pt(), weight);
        hComParton_[SL_topY]->Fill(partonTop1->p4().Rapidity(), weight);
        hComParton_[SL_topY]->Fill(partonTop2->p4().Rapidity(), weight);
        hComParton_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
        hComParton_[SL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()), weight);
        hComParton_[SL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()), weight);
        hComParton_[SL_ttbarPt]->Fill(partonTT.pt(), weight);
        hComParton_[SL_ttbarY]->Fill(partonTT.Rapidity(), weight);
        hComParton_[SL_ttbarMass]->Fill(partonTT.mass(), weight);
        hComParton_[SL_topPtTtbarSys]->Fill(partonTopPtAtCM, weight);

        hComPseudo_[SL_topPt]->Fill(pseudoTop1->pt(), weight);
        hComPseudo_[SL_topPt]->Fill(pseudoTop2->pt(), weight);
        hComPseudo_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), weight);
        hComPseudo_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), weight);
        hComPseudo_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), weight);
        hComPseudo_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), weight);
        hComPseudo_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), weight);
        hComPseudo_[SL_ttbarPt]->Fill(pseudoTT.pt(), weight);
        hComPseudo_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), weight);
        hComPseudo_[SL_ttbarMass]->Fill(pseudoTT.mass(), weight);
        hComPseudo_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, weight);

        // Fill response matrix no matter what parton level object acceptance is
        h2Com_[SL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt(), weight);
        h2Com_[SL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt(), weight);
        h2Com_[SL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
        h2Com_[SL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
        h2Com_[SL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
        h2Com_[SL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()), weight);
        h2Com_[SL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()), weight);
        h2Com_[SL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt(), weight);
        h2Com_[SL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity(), weight);
        h2Com_[SL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass(), weight);
        h2Com_[SL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM, weight);
      }
    }
  }
  else if ( (pseudoTopCh == 3 or pseudoTopCh == 4 or pseudoTopCh == 5 ) and
            isAcceptedFullLept(pseudoW11, pseudoW21, pseudoB1, pseudoB2) ) {
    h2PseudoChannel_->Fill(partonTopCh, pseudoTopCh, weight);

    hPseudo_[DL_topPt]->Fill(pseudoTop1->pt(), weight);
    hPseudo_[DL_topPt]->Fill(pseudoTop2->pt(), weight);
    hPseudo_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), weight);
    hPseudo_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), weight);
    hPseudo_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), weight);
    hPseudo_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), weight);
    hPseudo_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), weight);
    hPseudo_[DL_ttbarPt]->Fill(pseudoTT.pt(), weight);
    hPseudo_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), weight);
    hPseudo_[DL_ttbarMass]->Fill(pseudoTT.mass(), weight);
    hPseudo_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, weight);

    // Fill response matrix no matter what parton level object acceptance is
    h2_[DL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt(), weight);
    h2_[DL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt(), weight);
    h2_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
    h2_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
    h2_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
    h2_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()), weight);
    h2_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()), weight);
    h2_[DL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt(), weight);
    h2_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity(), weight);
    h2_[DL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass(), weight);
    h2_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM, weight);

    if ( channel == CH_FULLLEPTON ) {
      hChPseudo_Channel_->Fill(pseudoTopCh, weight);
      h2ChChannel_->Fill(partonTopCh, pseudoTopCh, weight);

      hChPseudo_[DL_topPt]->Fill(pseudoTop1->pt(), weight);
      hChPseudo_[DL_topPt]->Fill(pseudoTop2->pt(), weight);
      hChPseudo_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), weight);
      hChPseudo_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), weight);
      hChPseudo_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), weight);
      hChPseudo_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), weight);
      hChPseudo_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), weight);
      hChPseudo_[DL_ttbarPt]->Fill(pseudoTT.pt(), weight);
      hChPseudo_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), weight);
      hChPseudo_[DL_ttbarMass]->Fill(pseudoTT.mass(), weight);
      hChPseudo_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, weight);

      if ( isAcceptedFullLept(partonW11, partonW21, partonB1, partonB2) ) {
        hComParton_Channel_->Fill(partonTopCh, weight);
        hComPseudo_Channel_->Fill(pseudoTopCh, weight);
        h2ComChannel_->Fill(partonTopCh, pseudoTopCh, weight);

        hComParton_[DL_topPt]->Fill(partonTop1->pt(), weight);
        hComParton_[DL_topPt]->Fill(partonTop2->pt(), weight);
        hComParton_[DL_topY]->Fill(partonTop1->p4().Rapidity(), weight);
        hComParton_[DL_topY]->Fill(partonTop2->p4().Rapidity(), weight);
        hComParton_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
        hComParton_[DL_topPtLead]->Fill(std::max(partonTop1->pt(), partonTop2->pt()), weight);
        hComParton_[DL_topPtSubLead]->Fill(std::min(partonTop1->pt(), partonTop2->pt()), weight);
        hComParton_[DL_ttbarPt]->Fill(partonTT.pt(), weight);
        hComParton_[DL_ttbarY]->Fill(partonTT.Rapidity(), weight);
        hComParton_[DL_ttbarMass]->Fill(partonTT.mass(), weight);
        hComParton_[DL_topPtTtbarSys]->Fill(partonTopPtAtCM, weight);

        hComPseudo_[DL_topPt]->Fill(pseudoTop1->pt(), weight);
        hComPseudo_[DL_topPt]->Fill(pseudoTop2->pt(), weight);
        hComPseudo_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), weight);
        hComPseudo_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), weight);
        hComPseudo_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), weight);
        hComPseudo_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), weight);
        hComPseudo_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), weight);
        hComPseudo_[DL_ttbarPt]->Fill(pseudoTT.pt(), weight);
        hComPseudo_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), weight);
        hComPseudo_[DL_ttbarMass]->Fill(pseudoTT.mass(), weight);
        hComPseudo_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, weight);

        // Fill response matrix no matter what parton level object acceptance is
        h2Com_[DL_topPt]->Fill(pseudoTop1->pt(), partonTop1->pt(), weight);
        h2Com_[DL_topPt]->Fill(pseudoTop2->pt(), partonTop2->pt(), weight);
        h2Com_[DL_topY]->Fill(pseudoTop1->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
        h2Com_[DL_topY]->Fill(pseudoTop2->p4().Rapidity(), partonTop1->p4().Rapidity(), weight);
        h2Com_[DL_ttbarDelPhi]->Fill(reco::deltaPhi(pseudoTop1->phi(), pseudoTop2->phi()), reco::deltaPhi(partonTop1->phi(), partonTop2->phi()), weight);
        h2Com_[DL_topPtLead]->Fill(std::max(pseudoTop1->pt(), pseudoTop2->pt()), std::max(partonTop1->pt(), partonTop2->pt()), weight);
        h2Com_[DL_topPtSubLead]->Fill(std::min(pseudoTop1->pt(), pseudoTop2->pt()), std::min(partonTop1->pt(), partonTop2->pt()), weight);
        h2Com_[DL_ttbarPt]->Fill(pseudoTT.pt(), partonTT.pt(), weight);
        h2Com_[DL_ttbarY]->Fill(pseudoTT.Rapidity(), partonTT.Rapidity(), weight);
        h2Com_[DL_ttbarMass]->Fill(pseudoTT.mass(), partonTT.mass(), weight);
        h2Com_[DL_topPtTtbarSys]->Fill(pseudoTopPtAtCM, partonTopPtAtCM, weight);
      }
    }
  }

}

DEFINE_FWK_MODULE(CATGenTopAnalysis);

