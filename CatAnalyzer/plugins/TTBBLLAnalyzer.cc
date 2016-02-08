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

struct Histos
{
  typedef TH1D* H1;
  typedef TH2D* H2;

  H1 bjetsT_n, bjetsM_n, bjetsL_n;
  H1 jet1_btag, jet2_btag, jet3_btag, jet4_btag;

  H2 jet3_btag__jet4_btag;

  void book(TFileDirectory dir)
  {
    bjetsT_n = dir.make<TH1D>("bjetsT_n", "Tight b jet multiplicity;B jet multiplicity", 10, 0, 10);
    bjetsM_n = dir.make<TH1D>("bjetsM_n", "Medium b jet multiplicity;B jet multiplicity", 10, 0, 10);
    bjetsL_n = dir.make<TH1D>("bjetsL_n", "Loose b jet multiplicity;B jet multiplicity", 10, 0, 10);

    jet1_btag = dir.make<TH1D>("jet1_btag", "1st b discriminator;B discriminator", 10, 0, 1);
    jet2_btag = dir.make<TH1D>("jet2_btag", "2nd b discriminator;B discriminator", 10, 0, 1);
    jet3_btag = dir.make<TH1D>("jet3_btag", "3rd b discriminator;B discriminator", 10, 0, 1);
    jet4_btag = dir.make<TH1D>("jet4_btag", "4th b discriminator;B discriminator", 10, 0, 1);

    jet3_btag__jet4_btag = dir.make<TH2D>("jet3_btag__jet4_btag", "3rd vs 4th b discriminator;3rd b discrinimator;4th b discriminator", 10, 0, 1, 10, 0, 1);
  }
};

class TTBBLLAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  TTBBLLAnalyzer(const edm::ParameterSet& pset);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  const reco::Candidate* findTopMother(const reco::Candidate* p) const;
  typedef std::vector<int> vint;

  edm::EDGetTokenT<int> channelToken_;
  edm::EDGetTokenT<float> weightToken_;
  edm::EDGetTokenT<cat::LeptonCollection> leptonsToken_;
  edm::EDGetTokenT<cat::JetCollection> jetsToken_;
  edm::EDGetTokenT<float> metToken_, metphiToken_;

  //edm::EDGetTokenT<reco::GenParticleCollection> partonToken_;
  edm::EDGetTokenT<int> partonChannelToken_;
  edm::EDGetTokenT<vint> partonModesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> partonJetToken_;

  Histos heeS0_, heeS1_, heeS2_;
  Histos hmmS0_, hmmS1_, hmmS2_;
  Histos hemS0_, hemS1_, hemS2_;

  TTree* tree_;

  int b_channel_;
  float b_weight_;
  float b_met_pt_, b_met_phi_;

  float b_lep1_pt_, b_lep2_pt_, b_z_m_, b_z_pt_;

  float b_jet1_pt_, b_jet2_pt_, b_jet3_pt_, b_jet4_pt_;
  float b_jet1_eta_, b_jet2_eta_, b_jet3_eta_, b_jet4_eta_;
  float b_jet1_btag_, b_jet2_btag_, b_jet3_btag_, b_jet4_btag_;
  int b_jet1_hflav_, b_jet2_hflav_, b_jet3_hflav_, b_jet4_hflav_;
  int b_jet1_qflav_, b_jet2_qflav_, b_jet3_qflav_, b_jet4_qflav_;

  int b_bjetsT_n, b_bjetsM_n, b_bjetsL_n;
  float b_csvWeight_;

  int b_parton_channel_, b_parton_mode1_, b_parton_mode2_;
  int b_parton_jets20_n_, b_parton_jets30_n_;
  int b_parton_bjets20_n_, b_parton_bjets30_n_;

  float b_parton_l1_pt_, b_parton_l2_pt_;
  float b_parton_l1_eta_, b_parton_l2_eta_;
  float b_parton_topbjet1_pt_, b_parton_topbjet2_pt_;
  float b_parton_topbjet1_eta_, b_parton_topbjet2_eta_;
  int b_parton_addjets20_n_, b_parton_addjets30_n_, b_parton_addjets40_n_;
  int b_parton_addbjets20_n_, b_parton_addbjets30_n_, b_parton_addbjets40_n_;
  int b_parton_addcjets20_n_, b_parton_addcjets30_n_, b_parton_addcjets40_n_;

  const bool doTree_;
  const bool isTopMC_;

  CSVWeightEvaluator csvWeight_;
};

TTBBLLAnalyzer::TTBBLLAnalyzer(const edm::ParameterSet& pset):
  doTree_(pset.getParameter<bool>("doTree")),
  isTopMC_(pset.getParameter<bool>("isTopMC")),
  csvWeight_(CSVWeightEvaluator::ROOT)
{
  const auto srcLabel = pset.getParameter<edm::InputTag>("src");
  const auto srcLabelName = srcLabel.label();

  channelToken_ = consumes<int>(edm::InputTag(srcLabelName, "channel"));
  weightToken_ = consumes<float>(edm::InputTag(srcLabelName, "weight"));
  leptonsToken_ = consumes<cat::LeptonCollection>(edm::InputTag(srcLabelName, "leptons"));
  jetsToken_ = consumes<cat::JetCollection>(edm::InputTag(srcLabelName, "jets"));
  metToken_ = consumes<float>(edm::InputTag(srcLabelName, "met"));
  metphiToken_ = consumes<float>(edm::InputTag(srcLabelName, "metphi"));

  if ( isTopMC_ ) {
    const auto partonLabel = pset.getParameter<edm::InputTag>("partonTop");
    const auto partonLabelName = partonLabel.label();

    //partonToken_ = consumes<reco::GenParticleCollection>(partonLabel);
    partonChannelToken_ = consumes<int>(edm::InputTag(partonLabelName, "channel"));
    partonModesToken_ = consumes<vint>(edm::InputTag(partonLabelName, "modes"));
    partonJetToken_ = consumes<reco::GenJetCollection>(edm::InputTag(partonLabelName, "qcdJets"));
  }

  edm::Service<TFileService> fs;
  if ( doTree_ ) {
    tree_ = fs->make<TTree>("ttbb", "ttbb");

    tree_->Branch("channel", &b_channel_, "channel/I");
    tree_->Branch("weight", &b_weight_, "weight/F");

    tree_->Branch("met_pt", &b_met_pt_, "met_pt/F");
    tree_->Branch("met_phi", &b_met_phi_, "met_phi/F");

    tree_->Branch("lep1_pt", &b_lep1_pt_, "lep1_pt/F");
    tree_->Branch("lep2_pt", &b_lep2_pt_, "lep2_pt/F");
    tree_->Branch("z_m", &b_z_m_, "z_m/F");
    tree_->Branch("z_pt", &b_z_pt_, "z_pt/F");

    tree_->Branch("jet1_pt", &b_jet1_pt_, "jet1_pt/F");
    tree_->Branch("jet2_pt", &b_jet2_pt_, "jet2_pt/F");
    tree_->Branch("jet3_pt", &b_jet3_pt_, "jet3_pt/F");
    tree_->Branch("jet4_pt", &b_jet4_pt_, "jet4_pt/F");

    tree_->Branch("jet1_eta", &b_jet1_eta_, "jet1_eta/F");
    tree_->Branch("jet2_eta", &b_jet2_eta_, "jet2_eta/F");
    tree_->Branch("jet3_eta", &b_jet3_eta_, "jet3_eta/F");
    tree_->Branch("jet4_eta", &b_jet4_eta_, "jet4_eta/F");

    tree_->Branch("jet1_btag", &b_jet1_btag_, "jet1_btag/F");
    tree_->Branch("jet2_btag", &b_jet2_btag_, "jet2_btag/F");
    tree_->Branch("jet3_btag", &b_jet3_btag_, "jet3_btag/F");
    tree_->Branch("jet4_btag", &b_jet4_btag_, "jet4_btag/F");

    tree_->Branch("jet1_hflav", &b_jet1_hflav_, "jet1_hflav/I");
    tree_->Branch("jet2_hflav", &b_jet2_hflav_, "jet2_hflav/I");
    tree_->Branch("jet3_hflav", &b_jet3_hflav_, "jet3_hflav/I");
    tree_->Branch("jet4_hflav", &b_jet4_hflav_, "jet4_hflav/I");

    tree_->Branch("jet1_qflav", &b_jet1_qflav_, "jet1_qflav/I");
    tree_->Branch("jet2_qflav", &b_jet2_qflav_, "jet2_qflav/I");
    tree_->Branch("jet3_qflav", &b_jet3_qflav_, "jet3_qflav/I");
    tree_->Branch("jet4_qflav", &b_jet4_qflav_, "jet4_qflav/I");

    tree_->Branch("nbjetsT", &b_bjetsT_n, "nbjetsT/I");
    tree_->Branch("nbjetsM", &b_bjetsM_n, "nbjetsM/I");
    tree_->Branch("nbjetsL", &b_bjetsL_n, "nbjetsL/I");

    tree_->Branch("csvWeight", &b_csvWeight_, "csvWeight/F");

    if ( isTopMC_ ) {
      tree_->Branch("parton_channel", &b_parton_channel_, "parton_channel/I");
      tree_->Branch("parton_mode1", &b_parton_mode1_, "parton_mode1/I");
      tree_->Branch("parton_mode2", &b_parton_mode2_, "parton_mode2/I");
      tree_->Branch("parton_jets20_n", &b_parton_jets20_n_, "parton_jets20_n/I");
      tree_->Branch("parton_jets30_n", &b_parton_jets30_n_, "parton_jets30_n/I");
      tree_->Branch("parton_bjets20_n", &b_parton_bjets20_n_, "parton_bjets20_n/I");
      tree_->Branch("parton_bjets30_n", &b_parton_bjets30_n_, "parton_bjets30_n/I");

      tree_->Branch("parton_l1_pt", &b_parton_l1_pt_, "parton_l1_pt/F");
      tree_->Branch("parton_l2_pt", &b_parton_l2_pt_, "parton_l2_pt/F");
      tree_->Branch("parton_l1_eta", &b_parton_l1_eta_, "parton_l1_eta/F");
      tree_->Branch("parton_l2_eta", &b_parton_l2_eta_, "parton_l2_eta/F");

      tree_->Branch("parton_topbjet1_pt", &b_parton_topbjet1_pt_, "parton_topbjet1_pt/F");
      tree_->Branch("parton_topbjet2_pt", &b_parton_topbjet2_pt_, "parton_topbjet2_pt/F");
      tree_->Branch("parton_topbjet1_eta", &b_parton_topbjet1_eta_, "parton_topbjet1_eta/F");
      tree_->Branch("parton_topbjet2_eta", &b_parton_topbjet2_eta_, "parton_topbjet2_eta/F");

      tree_->Branch("parton_addjets20_n", &b_parton_addjets20_n_, "parton_addjets20_n/I");
      tree_->Branch("parton_addjets30_n", &b_parton_addjets30_n_, "parton_addjets30_n/I");
      tree_->Branch("parton_addjets40_n", &b_parton_addjets40_n_, "parton_addjets40_n/I");

      tree_->Branch("parton_addbjets20_n", &b_parton_addbjets20_n_, "parton_addbjets20_n/I");
      tree_->Branch("parton_addbjets30_n", &b_parton_addbjets30_n_, "parton_addbjets30_n/I");
      tree_->Branch("parton_addbjets40_n", &b_parton_addbjets40_n_, "parton_addbjets40_n/I");

      tree_->Branch("parton_addcjets20_n", &b_parton_addcjets20_n_, "parton_addcjets20_n/I");
      tree_->Branch("parton_addcjets30_n", &b_parton_addcjets30_n_, "parton_addcjets30_n/I");
      tree_->Branch("parton_addcjets40_n", &b_parton_addcjets40_n_, "parton_addcjets40_n/I");
    }
  }

  auto diree = fs->mkdir("ee");
  heeS0_.book(diree.mkdir("S0")); // Njet4
  heeS1_.book(diree.mkdir("S1")); // Medium B tag for leading two
  heeS2_.book(diree.mkdir("S2")); // Tight B tag for leading two

  auto dirmm = fs->mkdir("mm");
  hmmS0_.book(dirmm.mkdir("S0")); // Njet4
  hmmS1_.book(dirmm.mkdir("S1")); // Medium B tag for leading two
  hmmS2_.book(dirmm.mkdir("S2")); // Tight B tag for leading two

  auto direm = fs->mkdir("em");
  hemS0_.book(direm.mkdir("S0")); // Njet4
  hemS1_.book(direm.mkdir("S1")); // Medium B tag for leading two
  hemS2_.book(direm.mkdir("S2")); // Tight B tag for leading two
}

void TTBBLLAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&)
{
  // Initialize variables
  b_lep1_pt_ =  b_lep2_pt_ = b_z_m_ = b_z_pt_ = 0;

  b_jet1_pt_ = b_jet2_pt_ = b_jet3_pt_ = b_jet4_pt_ = 0;
  b_jet1_eta_ = b_jet2_eta_ = b_jet3_eta_ = b_jet4_eta_ = -999;
  b_jet1_btag_ = b_jet2_btag_ = b_jet3_btag_ = b_jet4_btag_ = -999;
  b_jet1_hflav_ = b_jet2_hflav_ = b_jet3_hflav_ = b_jet4_hflav_ = -999;
  b_jet1_qflav_ = b_jet2_qflav_ = b_jet3_qflav_ = b_jet4_qflav_ = -999;

  b_bjetsT_n = b_bjetsM_n = b_bjetsL_n = 0;

  if ( isTopMC_ ) {
    b_parton_channel_ = CH_NOTT;
    b_parton_mode1_ = b_parton_mode2_= CH_HADRON;
    b_parton_jets20_n_ = b_parton_jets30_n_= 0;
    b_parton_bjets20_n_ = b_parton_bjets30_n_= 0;

    b_parton_l1_pt_ = b_parton_l2_pt_ = 0;
    b_parton_l1_eta_ = b_parton_l2_eta_ = -999;

    b_parton_topbjet1_pt_ = b_parton_topbjet2_pt_ = 0;
    b_parton_topbjet1_eta_ = b_parton_topbjet2_eta_ = -999;

    b_parton_addjets20_n_ = b_parton_addjets30_n_ = b_parton_addjets40_n_ = 0;
    b_parton_addbjets20_n_ = b_parton_addbjets30_n_ = b_parton_addbjets40_n_ = 0;
    b_parton_addcjets20_n_ = b_parton_addcjets30_n_ = b_parton_addcjets40_n_ = 0;
  }

  // Start to read reco objects
  edm::Handle<int> iHandle;
  event.getByToken(channelToken_, iHandle);
  b_channel_ = *iHandle;

  edm::Handle<float> fHandle;
  event.getByToken(weightToken_, fHandle);
  b_weight_ = *fHandle;
  event.getByToken(metToken_, fHandle);
  b_met_pt_ = *fHandle;
  event.getByToken(metphiToken_, fHandle);
  b_met_phi_ = *fHandle;

  edm::Handle<cat::LeptonCollection> leptonsHandle;
  event.getByToken(leptonsToken_, leptonsHandle);
  if ( leptonsHandle->size() < 2 ) return; // Just for a confirmation
  b_lep1_pt_ = leptonsHandle->at(0).pt();
  b_lep2_pt_ = leptonsHandle->at(1).pt();
  const auto zP4 = leptonsHandle->at(0).p4()+leptonsHandle->at(1).p4();
  b_z_m_ = zP4.mass();
  b_z_pt_ = zP4.pt();

  edm::Handle<cat::JetCollection> jetsHandle;
  event.getByToken(jetsToken_, jetsHandle);
  const int nJets = jetsHandle->size();
  if ( nJets < 4 ) return;

  typedef std::pair<edm::Ref<cat::JetCollection>, double> JetRefQPair;
  std::vector<JetRefQPair> jetRefs;
  b_bjetsT_n = 0, b_bjetsM_n = 0, b_bjetsL_n = 0;
  for ( size_t i=0, n=jetsHandle->size(); i<n; ++i ) {
    const auto& jet = jetsHandle->at(i);
    const double bTag = jet.bDiscriminator(BTAG_CSVv2);
    if ( bTag > WP_BTAG_CSVv2L ) ++b_bjetsL_n;
    if ( bTag > WP_BTAG_CSVv2M ) ++b_bjetsM_n;
    if ( bTag > WP_BTAG_CSVv2T ) ++b_bjetsT_n;
    jetRefs.push_back(JetRefQPair(edm::Ref<cat::JetCollection>(jetsHandle, i), bTag));
  }
  std::sort(jetRefs.begin(), jetRefs.end(), [](const JetRefQPair& a, const JetRefQPair& b){return a.second > b.second;});
  b_csvWeight_  = csvWeight_(*jetRefs.at(0).first, CSVWeightEvaluator::CENTRAL);
  b_csvWeight_ *= csvWeight_(*jetRefs.at(1).first, CSVWeightEvaluator::CENTRAL);
  b_csvWeight_ *= csvWeight_(*jetRefs.at(2).first, CSVWeightEvaluator::CENTRAL);
  b_csvWeight_ *= csvWeight_(*jetRefs.at(3).first, CSVWeightEvaluator::CENTRAL);

  b_jet1_pt_ = jetRefs.at(0).first->pt();
  b_jet2_pt_ = jetRefs.at(1).first->pt();
  b_jet3_pt_ = jetRefs.at(2).first->pt();
  b_jet4_pt_ = jetRefs.at(3).first->pt();

  b_jet1_eta_ = jetRefs.at(0).first->eta();
  b_jet2_eta_ = jetRefs.at(1).first->eta();
  b_jet3_eta_ = jetRefs.at(2).first->eta();
  b_jet4_eta_ = jetRefs.at(3).first->eta();

  b_jet1_btag_ = jetRefs.at(0).second;
  b_jet2_btag_ = jetRefs.at(1).second;
  b_jet3_btag_ = jetRefs.at(2).second;
  b_jet4_btag_ = jetRefs.at(3).second;

  b_jet1_hflav_ = jetRefs.at(0).first->hadronFlavour();
  b_jet2_hflav_ = jetRefs.at(1).first->hadronFlavour();
  b_jet3_hflav_ = jetRefs.at(2).first->hadronFlavour();
  b_jet4_hflav_ = jetRefs.at(3).first->hadronFlavour();

  b_jet1_qflav_ = jetRefs.at(0).first->partonFlavour();
  b_jet2_qflav_ = jetRefs.at(1).first->partonFlavour();
  b_jet3_qflav_ = jetRefs.at(2).first->partonFlavour();
  b_jet4_qflav_ = jetRefs.at(3).first->partonFlavour();

  if ( isTopMC_ ) {
    // Move to the partons
    //edm::Handle<reco::GenParticleCollection> partonHandle;
    //event.getByToken(partonToken_, partonHandle);
    event.getByToken(partonChannelToken_, iHandle);
    b_parton_channel_ = *iHandle;
    edm::Handle<vint> viHandle;
    event.getByToken(partonModesToken_, viHandle);
    if ( viHandle->size() >= 1 ) b_parton_mode1_ = viHandle->at(0);
    if ( viHandle->size() >= 2 ) b_parton_mode2_ = viHandle->at(1);

    edm::Handle<reco::GenJetCollection> partonJetHandle;
    event.getByToken(partonJetToken_, partonJetHandle);
    std::vector<const reco::GenJet*> topbjets, addjets, addbjets, addcjets;
    for ( auto& jet : *partonJetHandle ) {
      if ( jet.pt() < 20 ) continue;
      if ( std::abs(jet.eta()) > 2.5 ) continue;

      bool isBjet = false, isCjet = false;
      bool isFromTop = false;
      for ( auto& con : jet.getGenConstituents() ) {
        if ( !isBjet and std::abs(con->pdgId()) == 5 ) isBjet = true;
        if ( !isCjet and std::abs(con->pdgId()) == 4 ) isCjet = true;
        if ( !isFromTop and findTopMother(con) ) isFromTop = true;
      }

      if ( isFromTop and isBjet ) topbjets.push_back(&jet);
      else if ( !isFromTop ) {
        addjets.push_back(&jet);
        if ( isBjet ) addbjets.push_back(&jet);
        if ( !isBjet and isCjet ) addcjets.push_back(&jet);
      }
    }

    if ( topbjets.size() >= 2 ) {
      std::nth_element(topbjets.begin(), topbjets.begin()+2, topbjets.end(),
                       [](const reco::GenJet* a, const reco::GenJet* b){return a->pt() > b->pt();});

      b_parton_topbjet1_pt_ = topbjets.at(0)->pt();
      b_parton_topbjet2_pt_ = topbjets.at(1)->pt();
      b_parton_topbjet1_eta_ = topbjets.at(0)->eta();
      b_parton_topbjet2_eta_ = topbjets.at(1)->eta();
    }

    for ( const auto& x : addjets ) {
      const double pt = x->pt();
      if ( pt > 20 ) ++b_parton_addjets20_n_;
      if ( pt > 30 ) ++b_parton_addjets30_n_;
      if ( pt > 40 ) ++b_parton_addjets40_n_;
    }

    for ( const auto& x : addbjets ) {
      const double pt = x->pt();
      if ( pt > 20 ) ++b_parton_addbjets20_n_;
      if ( pt > 30 ) ++b_parton_addbjets30_n_;
      if ( pt > 40 ) ++b_parton_addbjets40_n_;
    }

    for ( const auto& x : addcjets ) {
      const double pt = x->pt();
      if ( pt > 20 ) ++b_parton_addcjets20_n_;
      if ( pt > 30 ) ++b_parton_addcjets30_n_;
      if ( pt > 40 ) ++b_parton_addcjets40_n_;
    }
  }

  if ( b_channel_ == CH_MUMU ) { // CH_MUMU
    hmmS0_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
    hmmS0_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
    hmmS0_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

    hmmS0_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
    hmmS0_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
    hmmS0_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
    hmmS0_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

    hmmS0_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);

    if ( b_bjetsM_n >= 2 ) {
      hmmS1_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
      hmmS1_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
      hmmS1_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

      hmmS1_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hmmS1_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hmmS1_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hmmS1_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hmmS1_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }

    if ( b_bjetsT_n >= 2 ) {
      hmmS2_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
      hmmS2_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
      hmmS2_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

      hmmS2_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hmmS2_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hmmS2_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hmmS2_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hmmS2_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }
  }
  else if ( b_channel_ == CH_ELEL ) { // CH_ELEL
    heeS0_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
    heeS0_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
    heeS0_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

    heeS0_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
    heeS0_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
    heeS0_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
    heeS0_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

    heeS0_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);

    if ( b_bjetsM_n >= 2 ) {
      heeS1_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
      heeS1_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
      heeS1_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

      heeS1_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      heeS1_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      heeS1_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      heeS1_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      heeS1_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }

    if ( b_bjetsT_n >= 2 ) {
      heeS2_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
      heeS2_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
      heeS2_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

      heeS2_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      heeS2_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      heeS2_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      heeS2_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      heeS2_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }
  }
  else if ( b_channel_ == CH_MUEL ) { // CH_MUEL
    hemS0_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
    hemS0_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
    hemS0_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

    hemS0_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
    hemS0_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
    hemS0_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
    hemS0_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

    hemS0_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);

    if ( b_bjetsM_n >= 2 ) {
      hemS1_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
      hemS1_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
      hemS1_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

      hemS1_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hemS1_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hemS1_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hemS1_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hemS1_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }

    if ( b_bjetsT_n >= 2 ) {
      hemS2_.bjetsT_n->Fill(b_bjetsT_n, b_weight_);
      hemS2_.bjetsM_n->Fill(b_bjetsM_n, b_weight_);
      hemS2_.bjetsL_n->Fill(b_bjetsL_n, b_weight_);

      hemS2_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hemS2_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hemS2_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hemS2_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hemS2_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }
  }

  if ( doTree_ ) tree_->Fill();
}

const reco::Candidate* TTBBLLAnalyzer::findTopMother(const reco::Candidate* p) const
{
  if ( !p ) return 0;

  for ( size_t i=0, n=p->numberOfMothers(); i<n; ++i ) {
    const auto m = p->mother(i);
    if ( !m ) continue;
    if ( std::abs(m->pdgId()) == 6 ) return m;
    const auto mm = findTopMother(m);
    if ( mm ) return mm;
  }

  return 0;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTBBLLAnalyzer);

