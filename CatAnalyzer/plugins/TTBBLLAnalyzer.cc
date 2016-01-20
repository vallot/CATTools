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
    typedef std::vector<int> vint;

    edm::EDGetTokenT<int> channelToken_;
    edm::EDGetTokenT<vint> modesToken_;
    edm::EDGetTokenT<float> weightToken_;
    edm::EDGetTokenT<cat::LeptonCollection> leptonsToken_;
    edm::EDGetTokenT<cat::JetCollection> jetsToken_;
    edm::EDGetTokenT<float> metToken_, metphiToken_;

    Histos heeS0_, heeS1_, heeS2_;
    Histos hmmS0_, hmmS1_, hmmS2_;
    Histos hemS0_, hemS1_, hemS2_;

    TTree* tree_;
    int b_channel_, b_mode1_, b_mode2_;
    float b_weight_;
    float b_met_pt_, b_met_phi_;

    float b_lep1_pt_, b_lep2_pt_, b_z_m_, b_z_pt_;

    float b_jet1_pt_, b_jet2_pt_, b_jet3_pt_, b_jet4_pt_;
    float b_jet1_btag_, b_jet2_btag_, b_jet3_btag_, b_jet4_btag_;
    int b_jet1_hflav_, b_jet2_hflav_, b_jet3_hflav_, b_jet4_hflav_;
    int b_jet1_qflav_, b_jet2_qflav_, b_jet3_qflav_, b_jet4_qflav_;

    int b_nbjetsT_, b_nbjetsM_, b_nbjetsL_;
};

TTBBLLAnalyzer::TTBBLLAnalyzer(const edm::ParameterSet& pset)
{
  const auto srcLabel = pset.getParameter<edm::InputTag>("src");
  const auto srcLabelName = srcLabel.label();

  channelToken_ = consumes<int>(edm::InputTag(srcLabelName, "channel"));
  modesToken_ = consumes<vint>(edm::InputTag(srcLabelName, "modes"));
  weightToken_ = consumes<float>(edm::InputTag(srcLabelName, "weight"));
  leptonsToken_ = consumes<cat::LeptonCollection>(edm::InputTag(srcLabelName, "leptons"));
  jetsToken_ = consumes<cat::JetCollection>(edm::InputTag(srcLabelName, "jets"));
  metToken_ = consumes<float>(edm::InputTag(srcLabelName, "met"));
  metphiToken_ = consumes<float>(edm::InputTag(srcLabelName, "metphi"));

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ttbb", "ttbb");

  tree_->Branch("channel", &b_channel_, "channel/I");
  tree_->Branch("mode1", &b_mode1_, "mode1/I");
  tree_->Branch("mode2", &b_mode2_, "mode2/I");
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

  tree_->Branch("nbjetsT", &b_nbjetsT_, "nbjetsT/I");
  tree_->Branch("nbjetsM", &b_nbjetsM_, "nbjetsM/I");
  tree_->Branch("nbjetsL", &b_nbjetsL_, "nbjetsL/I");

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
  edm::Handle<int> iHandle;
  event.getByToken(channelToken_, iHandle);
  b_channel_ = *iHandle;
  edm::Handle<vint> viHandle;
  event.getByToken(modesToken_, viHandle);
  b_mode1_ = b_mode2_ = -1;
  if ( viHandle->size() >= 1 ) b_mode1_ = viHandle->at(0);
  if ( viHandle->size() >= 2 ) b_mode2_ = viHandle->at(1);

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
  b_nbjetsT_ = 0, b_nbjetsM_ = 0, b_nbjetsL_ = 0;
  for ( size_t i=0, n=jetsHandle->size(); i<n; ++i ) {
    const auto& jet = jetsHandle->at(i);
    const double bTag = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    if ( bTag > 0.605 ) ++b_nbjetsL_;
    if ( bTag > 0.890 ) ++b_nbjetsM_;
    if ( bTag > 0.970 ) ++b_nbjetsT_;
    jetRefs.push_back(JetRefQPair(edm::Ref<cat::JetCollection>(jetsHandle, i), bTag));
  }
  std::sort(jetRefs.begin(), jetRefs.end(), [](const JetRefQPair& a, const JetRefQPair& b){return a.second > b.second;});
  b_jet1_pt_ = jetRefs.at(0).first->pt();
  b_jet2_pt_ = jetRefs.at(1).first->pt();
  b_jet3_pt_ = jetRefs.at(2).first->pt();
  b_jet4_pt_ = jetRefs.at(3).first->pt();

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

  if ( b_channel_ == 0 ) { // CH_MUMU
    hmmS0_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
    hmmS0_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
    hmmS0_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

    hmmS0_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
    hmmS0_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
    hmmS0_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
    hmmS0_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

    hmmS0_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);

    if ( b_nbjetsM_ >= 2 ) {
      hmmS1_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
      hmmS1_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
      hmmS1_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

      hmmS1_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hmmS1_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hmmS1_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hmmS1_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hmmS1_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }

    if ( b_nbjetsT_ >= 2 ) {
      hmmS2_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
      hmmS2_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
      hmmS2_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

      hmmS2_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hmmS2_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hmmS2_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hmmS2_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hmmS2_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }
  }
  else if ( b_channel_ == 1 ) { // CH_ELEL
    heeS0_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
    heeS0_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
    heeS0_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

    heeS0_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
    heeS0_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
    heeS0_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
    heeS0_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

    heeS0_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);

    if ( b_nbjetsM_ >= 2 ) {
      heeS1_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
      heeS1_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
      heeS1_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

      heeS1_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      heeS1_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      heeS1_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      heeS1_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      heeS1_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }

    if ( b_nbjetsT_ >= 2 ) {
      heeS2_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
      heeS2_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
      heeS2_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

      heeS2_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      heeS2_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      heeS2_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      heeS2_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      heeS2_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }
  }
  else if ( b_channel_ == 2 ) { // CH_MUEL
    hemS0_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
    hemS0_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
    hemS0_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

    hemS0_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
    hemS0_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
    hemS0_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
    hemS0_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

    hemS0_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);

    if ( b_nbjetsM_ >= 2 ) {
      hemS1_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
      hemS1_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
      hemS1_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

      hemS1_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hemS1_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hemS1_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hemS1_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hemS1_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }

    if ( b_nbjetsT_ >= 2 ) {
      hemS2_.bjetsT_n->Fill(b_nbjetsT_, b_weight_);
      hemS2_.bjetsM_n->Fill(b_nbjetsM_, b_weight_);
      hemS2_.bjetsL_n->Fill(b_nbjetsL_, b_weight_);

      hemS2_.jet1_btag->Fill(b_jet1_btag_, b_weight_);
      hemS2_.jet2_btag->Fill(b_jet2_btag_, b_weight_);
      hemS2_.jet3_btag->Fill(b_jet3_btag_, b_weight_);
      hemS2_.jet4_btag->Fill(b_jet4_btag_, b_weight_);

      hemS2_.jet3_btag__jet4_btag->Fill(b_jet3_btag_, b_jet4_btag_, b_weight_);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTBBLLAnalyzer);

