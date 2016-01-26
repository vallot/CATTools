#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CatAnalyzer/interface/BTagScaleFactorEvaluators.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

using namespace std;
using namespace cat;

class TTLLAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  TTLLAnalyzer(const edm::ParameterSet& pset);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // Objects from the TTLLEventSelector
  edm::EDGetTokenT<int> channelToken_;
  edm::EDGetTokenT<float> weightToken_;
  edm::EDGetTokenT<cat::LeptonCollection> leptonsToken_;
  edm::EDGetTokenT<cat::JetCollection> jetsToken_;
  edm::EDGetTokenT<float> metToken_;

  // Objects from the TTLLKinSolutionProducer
  typedef std::vector<float> vfloat;
  edm::EDGetTokenT<edm::View<reco::Candidate> > candsToken_;
  edm::EDGetTokenT<vfloat> auxToken_;

  // Objects from the Parton level
  edm::EDGetTokenT<reco::GenParticleCollection> partonTopToken_;
  edm::EDGetTokenT<int> partonTopChToken_;

  // Objects from the pseudo top
  edm::EDGetTokenT<reco::GenParticleCollection> pseudoTopToken_;

  typedef TH1D* H1;
  typedef TH2D* H2;

  TTree* tree_;

  int b_channel_;
  float b_weight_;
  float b_met_pt_;

  float b_j1_pt_, b_j2_pt_, b_j3_pt_, b_j4_pt_;
  float b_j1_btag_, b_j2_btag_, b_j3_btag_, b_j4_btag_;
  float b_csvWeight_;

  float b_kinfit_quality_;
  float b_l1_pt_, b_l2_pt_;
  float b_b1_pt_, b_b2_pt_;
  float b_t1_pt_, b_t1_y_, b_t1_m_;
  float b_t2_pt_, b_t2_y_, b_t2_m_;
  float b_tt_pt_, b_tt_y_, b_tt_m_;
  float b_tt_dphi_;

  int b_partonTopChannel_;
  float b_partonL1_pt_, b_partonL2_pt_;
  float b_partonB1_pt_, b_partonB2_pt_;
  float b_partonT1_pt_, b_partonT1_y_, b_partonT1_m_;
  float b_partonT2_pt_, b_partonT2_y_, b_partonT2_m_;
  float b_partonTT_pt_, b_partonTT_y_, b_partonTT_m_;
  float b_partonTT_dphi_;

  int b_pseudoTopChannel_;
  float b_pseudoL1_pt_, b_pseudoL2_pt_;
  float b_pseudoB1_pt_, b_pseudoB2_pt_;
  float b_pseudoT1_pt_, b_pseudoT1_y_, b_pseudoT1_m_;
  float b_pseudoT2_pt_, b_pseudoT2_y_, b_pseudoT2_m_;
  float b_pseudoTT_pt_, b_pseudoTT_y_, b_pseudoTT_m_;
  float b_pseudoTT_dphi_;

  const bool doTree_;
  const bool isTopMC_;
  const std::string bTagName_;

  CSVWeightEvaluator csvWeight_;
};

TTLLAnalyzer::TTLLAnalyzer(const edm::ParameterSet& pset):
  doTree_(pset.getParameter<bool>("doTree")),
  isTopMC_(pset.getParameter<bool>("isTopMC")),
  bTagName_(pset.getParameter<std::string>("bTagName"))
{
  const auto recoLabel = pset.getParameter<edm::InputTag>("recoObjects").label();
  channelToken_ = consumes<int>(edm::InputTag(recoLabel, "channel"));
  weightToken_ = consumes<float>(edm::InputTag(recoLabel, "weight"));
  leptonsToken_ = consumes<cat::LeptonCollection>(edm::InputTag(recoLabel, "leptons"));
  jetsToken_ = consumes<cat::JetCollection>(edm::InputTag(recoLabel, "jets"));
  metToken_ = consumes<float>(edm::InputTag(recoLabel, "met"));

  const auto kinLabel = pset.getParameter<edm::InputTag>("kinfit");
  candsToken_ = consumes<edm::View<reco::Candidate> >(kinLabel);
  auxToken_ = consumes<vfloat>(edm::InputTag(kinLabel.label(), "aux"));

  if ( isTopMC_ ) {
    const auto partonTopLabel = pset.getParameter<edm::InputTag>("partonTop");
    partonTopToken_ = consumes<reco::GenParticleCollection>(partonTopLabel);
    partonTopChToken_ = consumes<int>(edm::InputTag(partonTopLabel.label(), "channel"));

    const auto pseudoTopLabel = pset.getParameter<edm::InputTag>("pseudoTop");
    pseudoTopToken_ = consumes<reco::GenParticleCollection>(pseudoTopLabel);
  }

  edm::Service<TFileService> fs;
  if ( doTree_ ) {
    tree_ = fs->make<TTree>("event", "event");
    tree_->Branch("weight", &b_weight_, "weight/F");
    tree_->Branch("channel", &b_channel_, "channel/I");
    tree_->Branch("met_pt", &b_met_pt_, "met_pt/F");

    tree_->Branch("j1_pt", &b_j1_pt_, "j1_pt/F");
    tree_->Branch("j2_pt", &b_j2_pt_, "j2_pt/F");
    tree_->Branch("j3_pt", &b_j3_pt_, "j3_pt/F");
    tree_->Branch("j4_pt", &b_j4_pt_, "j4_pt/F");

    tree_->Branch("j1_btag", &b_j1_btag_, "j1_btag/F");
    tree_->Branch("j2_btag", &b_j2_btag_, "j2_btag/F");
    tree_->Branch("j3_btag", &b_j3_btag_, "j3_btag/F");
    tree_->Branch("j4_btag", &b_j4_btag_, "j4_btag/F");

    tree_->Branch("csvWeight", &b_csvWeight_, "csvWeight/F");

    tree_->Branch("kinfit_quality", &b_kinfit_quality_, "kinfit_quality/F");
    tree_->Branch("l1_pt", &b_l1_pt_, "l1_pt/F");
    tree_->Branch("l2_pt", &b_l2_pt_, "l2_pt/F");
    tree_->Branch("b1_pt", &b_b1_pt_, "b1_pt/F");
    tree_->Branch("b2_pt", &b_b2_pt_, "b2_pt/F");
    tree_->Branch("t1_pt", &b_t1_pt_, "t1_pt/F");
    tree_->Branch("t2_pt", &b_t2_pt_, "t2_pt/F");
    tree_->Branch("t1_y", &b_t1_y_, "t1_y/F");
    tree_->Branch("t2_y", &b_t2_y_, "t2_y/F");
    tree_->Branch("t1_m", &b_t1_m_, "t1_m/F");
    tree_->Branch("t2_m", &b_t2_m_, "t2_m/F");
    tree_->Branch("tt_pt", &b_tt_pt_, "tt_pt/F");
    tree_->Branch("tt_y", &b_tt_y_, "tt_y/F");
    tree_->Branch("tt_m", &b_tt_m_, "tt_m/F");
    tree_->Branch("tt_dphi", &b_tt_dphi_, "tt_dphi/F");

    if ( isTopMC_ ) {
      tree_->Branch("partonTopChannel", &b_partonTopChannel_, "partonTopChannel/I");
      tree_->Branch("partonL1_pt"  , &b_partonL1_pt_  , "partonL1_pt/F"  );
      tree_->Branch("partonL2_pt"  , &b_partonL2_pt_  , "partonL2_pt/F"  );
      tree_->Branch("partonB1_pt"  , &b_partonB1_pt_  , "partonB1_pt/F"  );
      tree_->Branch("partonB2_pt"  , &b_partonB2_pt_  , "partonB2_pt/F"  );
      tree_->Branch("partonT1_pt"  , &b_partonT1_pt_  , "partonT1_pt/F"  );
      tree_->Branch("partonT2_pt"  , &b_partonT2_pt_  , "partonT2_pt/F"  );
      tree_->Branch("partonT1_y"   , &b_partonT1_y_   , "partonT1_y/F"   );
      tree_->Branch("partonT2_y"   , &b_partonT2_y_   , "partonT2_y/F"   );
      tree_->Branch("partonT1_m"   , &b_partonT1_m_   , "partonT1_m/F"   );
      tree_->Branch("partonT2_m"   , &b_partonT2_m_   , "partonT2_m/F"   );
      tree_->Branch("partonTT_pt"  , &b_partonTT_pt_  , "partonTT_pt/F"  );
      tree_->Branch("partonTT_y"   , &b_partonTT_y_   , "partonTT_y/F"   );
      tree_->Branch("partonTT_m"   , &b_partonTT_m_   , "partonTT_m/F"   );
      tree_->Branch("partonTT_dphi", &b_partonTT_dphi_, "partonTT_dphi/F");

      tree_->Branch("pseudoTopChannel", &b_pseudoTopChannel_, "pseudoTopChannel/I");
      tree_->Branch("pseudoL1_pt"  , &b_pseudoL1_pt_  , "pseudoL1_pt/F"  );
      tree_->Branch("pseudoL2_pt"  , &b_pseudoL2_pt_  , "pseudoL2_pt/F"  );
      tree_->Branch("pseudoB1_pt"  , &b_pseudoB1_pt_  , "pseudoB1_pt/F"  );
      tree_->Branch("pseudoB2_pt"  , &b_pseudoB2_pt_  , "pseudoB2_pt/F"  );
      tree_->Branch("pseudoT1_pt"  , &b_pseudoT1_pt_  , "pseudoT1_pt/F"  );
      tree_->Branch("pseudoT2_pt"  , &b_pseudoT2_pt_  , "pseudoT2_pt/F"  );
      tree_->Branch("pseudoT1_y"   , &b_pseudoT1_y_   , "pseudoT1_y/F"   );
      tree_->Branch("pseudoT2_y"   , &b_pseudoT2_y_   , "pseudoT2_y/F"   );
      tree_->Branch("pseudoT1_m"   , &b_pseudoT1_m_   , "pseudoT1_m/F"   );
      tree_->Branch("pseudoT2_m"   , &b_pseudoT2_m_   , "pseudoT2_m/F"   );
      tree_->Branch("pseudoTT_pt"  , &b_pseudoTT_pt_  , "pseudoTT_pt/F"  );
      tree_->Branch("pseudoTT_y"   , &b_pseudoTT_y_   , "pseudoTT_y/F"   );
      tree_->Branch("pseudoTT_m"   , &b_pseudoTT_m_   , "pseudoTT_m/F"   );
      tree_->Branch("pseudoTT_dphi", &b_pseudoTT_dphi_, "pseudoTT_dphi/F");
    }
  }
}

void TTLLAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&)
{
  // Initialize variables
  b_j1_pt_ = b_j2_pt_ = b_j3_pt_ = b_j4_pt_ = 0;
  b_j1_btag_ = b_j2_btag_ = b_j3_btag_ = b_j4_btag_ = -999;
  b_csvWeight_ = 1.;

  b_kinfit_quality_ = -1;
  b_l1_pt_ = b_l2_pt_ = b_b1_pt_ = b_b2_pt_ = 0;
  b_t1_pt_ = 0; b_t1_y_ = -999; b_t1_m_ = 0;
  b_t2_pt_ = 0; b_t2_y_ = -999; b_t2_m_ = 0;
  b_tt_pt_ = 0; b_tt_y_ = -999; b_tt_m_ = 0;
  b_tt_dphi_ = -999;

  b_partonTopChannel_ = CH_NOTT;
  b_partonL1_pt_ = b_partonL2_pt_ = 0;
  b_partonB1_pt_ = b_partonB2_pt_ = 0;
  b_partonT1_pt_ = 0; b_partonT1_y_ = -999; b_partonT1_m_ = 0;
  b_partonT2_pt_ = 0; b_partonT2_y_ = -999; b_partonT2_m_ = 0;
  b_partonTT_pt_ = 0; b_partonTT_y_ = -999; b_partonTT_m_ = 0;
  b_partonTT_dphi_ = -999;

  b_pseudoTopChannel_ = CH_NOTT;
  b_pseudoL1_pt_ = b_pseudoL2_pt_ = 0;
  b_pseudoB1_pt_ = b_pseudoB2_pt_ = 0;
  b_pseudoT1_pt_ = 0; b_pseudoT1_y_ = -999; b_pseudoT1_m_ = 0;
  b_pseudoT2_pt_ = 0; b_pseudoT2_y_ = -999; b_pseudoT2_m_ = 0;
  b_pseudoTT_pt_ = 0; b_pseudoTT_y_ = -999; b_pseudoTT_m_ = 0;
  b_pseudoTT_dphi_ = -999;

  edm::Handle<int> iHandle;
  event.getByToken(channelToken_, iHandle);
  b_channel_ = *iHandle;
  if ( b_channel_ != CH_FULLLEPTON ) return; // Just for a confirmation

  edm::Handle<float> fHandle;
  event.getByToken(weightToken_, fHandle);
  b_weight_ = *fHandle;
  event.getByToken(metToken_, fHandle);
  b_met_pt_ = *fHandle;

  edm::Handle<cat::LeptonCollection> leptonsHandle;
  event.getByToken(leptonsToken_, leptonsHandle);
  if ( leptonsHandle->size() < 2 ) return;
  const auto& lepton1 = leptonsHandle->at(0);
  const auto& lepton2 = leptonsHandle->at(1);
  switch ( std::abs(lepton1.pdgId())+std::abs(lepton2.pdgId()) ) {
    case 11+11: b_channel_ = CH_ELEL; break;
    case 11+13: b_channel_ = CH_MUEL; break;
    case 13+13: b_channel_ = CH_MUMU; break;
    default: b_channel_ = CH_NOLL;
  }
  b_l1_pt_ = lepton1.pt();
  b_l2_pt_ = lepton2.pt();

  edm::Handle<cat::JetCollection> jetsHandle;
  event.getByToken(jetsToken_, jetsHandle);
  const int nJet = jetsHandle->size();
  if ( nJet > 0 ) {
    b_j1_pt_ = jetsHandle->at(0).pt();
    b_j1_btag_ = jetsHandle->at(0).bDiscriminator(bTagName_);
  }
  if ( nJet > 1 ) {
    b_j2_pt_ = jetsHandle->at(1).pt();
    b_j2_btag_ = jetsHandle->at(1).bDiscriminator(bTagName_);
  }
  if ( nJet > 2 ) {
    b_j3_pt_ = jetsHandle->at(2).pt();
    b_j3_btag_ = jetsHandle->at(2).bDiscriminator(bTagName_);
  }
  if ( nJet > 3 ) {
    b_j4_pt_ = jetsHandle->at(3).pt();
    b_j4_btag_ = jetsHandle->at(3).bDiscriminator(bTagName_);
  }
  for ( const auto& jet : *jetsHandle ) b_csvWeight_ *= csvWeight_(jet, CSVWeightEvaluator::CENTRAL);

  edm::Handle<edm::View<reco::Candidate> > candsHandle;
  event.getByToken(candsToken_, candsHandle);
  if ( candsHandle->size() >= 7 ) {
    edm::Handle<vfloat> vfHandle;
    event.getByToken(auxToken_, vfHandle);
    b_kinfit_quality_ = vfHandle->at(0);

    const auto& tt = candsHandle->at(0);
    const auto& t1 = candsHandle->at(1);
    const auto& t2 = candsHandle->at(2);
    const auto& w1 = candsHandle->at(3);
    const auto& w2 = candsHandle->at(4);
    //const auto& n1 = candsHandle->at(5);
    //const auto& n2 = cnadsHandle->at(6);

    const auto b1 = t1.p4()-w1.p4();
    const auto b2 = t2.p4()-w2.p4();
    b_b1_pt_ = b1.pt();
    b_b2_pt_ = b2.pt();
    b_t1_pt_ = t1.pt(); b_t1_y_ = t1.p4().Rapidity(); b_t1_m_ = t1.mass();
    b_t2_pt_ = t2.pt(); b_t2_y_ = t2.p4().Rapidity(); b_t2_m_ = t2.mass();
    b_tt_pt_ = tt.pt(); b_tt_y_ = tt.p4().Rapidity(); b_tt_m_ = tt.mass();
    b_tt_dphi_ = std::abs(deltaPhi(t1.phi(), t2.phi()));
  }

  if ( isTopMC_ ) {
    edm::Handle<reco::GenParticleCollection> partonTopHandle;
    event.getByToken(partonTopToken_, partonTopHandle);
    event.getByToken(partonTopChToken_, iHandle);
    b_partonTopChannel_ = *iHandle;

    edm::Handle<reco::GenParticleCollection> pseudoTopHandle;
    event.getByToken(pseudoTopToken_, pseudoTopHandle);
    b_pseudoTopChannel_ = CH_NOTT;

    auto gtByPtPtr = [](const reco::GenParticle* a, const reco::GenParticle* b) { return a->pt() > b->pt(); };

    do {
      if ( partonTopHandle->empty() ) break;
      std::vector<const reco::GenParticle*> topquarks;
      for ( auto& x : *partonTopHandle ) {
        if ( std::abs(x.pdgId()) == 6 ) topquarks.push_back(&x);
      }
      if ( topquarks.size() < 2 ) break;
      std::nth_element(topquarks.begin(), topquarks.begin()+2, topquarks.end(), gtByPtPtr);

      const auto t1 = topquarks.at(0);
      const auto t2 = topquarks.at(1);
      if ( t1->numberOfDaughters() < 2 or t2->numberOfDaughters() < 2 ) break;

      const auto w1 = t1->daughter(0);
      const auto w2 = t2->daughter(0);
      if ( !w1 or !w2 ) break;

      const auto b1 = t1->daughter(1);
      const auto b2 = t2->daughter(1);
      if ( !b1 or !b2 ) break;

      auto l1 = w1->daughter(0);
      auto l2 = w2->daughter(0);
      if ( !l1 or !l2 ) break;
      if ( l1->numberOfDaughters() > 1 ) l1 = l1->daughter(0);
      if ( l2->numberOfDaughters() > 1 ) l2 = l2->daughter(0);

      const int sumId = std::abs(l1->pdgId()) + std::abs(l2->pdgId());
      if      ( sumId < 11 ) b_partonTopChannel_ = CH_FULLHADRON;
      else if ( sumId < 15 ) b_partonTopChannel_ = CH_SEMILEPTON;

      if ( sumId != 11+11 and sumId != 11+13 and sumId != 13+13 ) break;

      b_partonT1_pt_ = t1->pt();
      b_partonT2_pt_ = t2->pt();
      b_partonT1_y_ = t1->p4().Rapidity();
      b_partonT2_y_ = t2->p4().Rapidity();
      b_partonT1_m_ = t1->mass();
      b_partonT2_m_ = t2->mass();

      const auto ttP4 = t1->p4()+t2->p4();
      b_partonTT_pt_ = ttP4.pt();
      b_partonTT_y_ = ttP4.Rapidity();
      b_partonTT_m_ = ttP4.mass();
      b_partonTT_dphi_ = std::abs(deltaPhi(t1->phi(), t2->phi()));

      b_partonL1_pt_ = l1->pt();
      b_partonL2_pt_ = l2->pt();
      b_partonB1_pt_ = b1->pt();
      b_partonB2_pt_ = b2->pt();
    } while ( false );

    do {
      if ( pseudoTopHandle->empty() ) break;
      std::vector<const reco::GenParticle*> topquarks;
      for ( auto& x : *pseudoTopHandle ) {
        if ( std::abs(x.pdgId()) == 6 ) topquarks.push_back(&x);
      }
      if ( topquarks.size() < 2 ) break;
      std::nth_element(topquarks.begin(), topquarks.begin()+2, topquarks.end(), gtByPtPtr);

      const auto t1 = topquarks.at(0);
      const auto t2 = topquarks.at(1);
      if ( t1->numberOfDaughters() < 2 or t2->numberOfDaughters() < 2 ) break;

      const auto w1 = t1->daughter(0);
      const auto w2 = t2->daughter(0);
      if ( !w1 or !w2 ) break;

      const auto b1 = t1->daughter(1);
      const auto b2 = t2->daughter(1);
      if ( !b1 or !b2 ) break;

      const auto l1 = w1->daughter(0);
      const auto l2 = w2->daughter(0);
      if ( !l1 or !l2 ) break;

      const int sumId = std::abs(l1->pdgId()) + std::abs(l2->pdgId());
      if      ( sumId < 11 ) b_pseudoTopChannel_ = CH_FULLHADRON;
      else if ( sumId < 15 ) b_pseudoTopChannel_ = CH_SEMILEPTON;

      if ( sumId != 11+11 and sumId != 11+13 and sumId != 13+13 ) break;
      b_pseudoTopChannel_ = CH_FULLLEPTON;

      b_pseudoT1_pt_ = t1->pt();
      b_pseudoT2_pt_ = t2->pt();
      b_pseudoT1_y_ = t1->p4().Rapidity();
      b_pseudoT2_y_ = t2->p4().Rapidity();
      b_pseudoT1_m_ = t1->mass();
      b_pseudoT2_m_ = t2->mass();

      const auto ttP4 = t1->p4()+t2->p4();
      b_pseudoTT_pt_ = ttP4.pt();
      b_pseudoTT_y_ = ttP4.Rapidity();
      b_pseudoTT_m_ = ttP4.mass();
      b_pseudoTT_dphi_ = std::abs(deltaPhi(t1->phi(), t2->phi()));

      b_pseudoL1_pt_ = l1->pt();
      b_pseudoL2_pt_ = l2->pt();
      b_pseudoB1_pt_ = b1->pt();
      b_pseudoB2_pt_ = b2->pt();
    } while ( false );
  }

  if ( doTree_ ) tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLLAnalyzer);

