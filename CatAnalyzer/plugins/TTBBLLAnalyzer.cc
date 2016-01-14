#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "TH1D.h"
#include "TH2D.h"

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
    edm::EDGetTokenT<int> channelToken_;
    edm::EDGetTokenT<float> weightToken_;
    edm::EDGetTokenT<cat::LeptonCollection> leptonsToken_;
    edm::EDGetTokenT<cat::JetCollection> jetsToken_;
    edm::EDGetTokenT<float> metToken_, metphiToken_;

    Histos heeS0_, heeS1_, heeS2_;
    Histos hmmS0_, hmmS1_, hmmS2_;
    Histos hemS0_, hemS1_, hemS2_;

};

TTBBLLAnalyzer::TTBBLLAnalyzer(const edm::ParameterSet& pset)
{
  const auto srcLabel = pset.getParameter<edm::InputTag>("src");
  const auto srcLabelName = srcLabel.label();

  channelToken_ = consumes<int>(edm::InputTag(srcLabelName, "channel"));
  weightToken_ = consumes<float>(edm::InputTag(srcLabelName, "weight"));
  leptonsToken_ = consumes<cat::LeptonCollection>(edm::InputTag(srcLabelName, "leptons"));
  jetsToken_ = consumes<cat::JetCollection>(edm::InputTag(srcLabelName, "jets"));
  metToken_ = consumes<float>(edm::InputTag(srcLabelName, "met"));
  metphiToken_ = consumes<float>(edm::InputTag(srcLabelName, "metphi"));

  edm::Service<TFileService> fs;
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
  const int channel = *iHandle;

  edm::Handle<float> fHandle;
  event.getByToken(weightToken_, fHandle);
  const double weight = *fHandle;
  //event.getByToken(metToken_, fHandle);
  //const double met_pt = *fHandle;
  //event.getByToken(metphiToken_, fHandle);
  //const double met_phi = *fHandle;

  //edm::Handle<cat::LeptonCollection> leptonsHandle;
  //event.getByToken(leptonsToken_, leptonsHandle);

  edm::Handle<cat::JetCollection> jetsHandle;
  event.getByToken(jetsToken_, jetsHandle);
  const int nJets = jetsHandle->size();
  if ( nJets < 4 ) return;

  std::vector<double> bTags;
  int nBjetsT = 0, nBjetsM = 0, nBjetsL = 0;
  for ( auto& jet : *jetsHandle ) {
    const double bTag = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    if ( bTag > 0.605 ) ++nBjetsL;
    if ( bTag > 0.890 ) ++nBjetsM;
    if ( bTag > 0.970 ) ++nBjetsT;
    bTags.push_back(bTag);
  }
  std::sort(bTags.begin(), bTags.end(), [](const double a, const double b){return a > b;});

  if ( channel == 0 ) { // CH_MUMU
    hmmS0_.bjetsT_n->Fill(nBjetsT, weight);
    hmmS0_.bjetsM_n->Fill(nBjetsM, weight);
    hmmS0_.bjetsL_n->Fill(nBjetsL, weight);

    hmmS0_.jet1_btag->Fill(bTags[0], weight);
    hmmS0_.jet2_btag->Fill(bTags[1], weight);
    hmmS0_.jet3_btag->Fill(bTags[2], weight);
    hmmS0_.jet4_btag->Fill(bTags[3], weight);

    hmmS0_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);

    if ( nBjetsM >= 2 ) {
      hmmS1_.bjetsT_n->Fill(nBjetsT, weight);
      hmmS1_.bjetsM_n->Fill(nBjetsM, weight);
      hmmS1_.bjetsL_n->Fill(nBjetsL, weight);

      hmmS1_.jet1_btag->Fill(bTags[0], weight);
      hmmS1_.jet2_btag->Fill(bTags[1], weight);
      hmmS1_.jet3_btag->Fill(bTags[2], weight);
      hmmS1_.jet4_btag->Fill(bTags[3], weight);

      hmmS1_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);
    }

    if ( nBjetsT >= 2 ) {
      hmmS2_.bjetsT_n->Fill(nBjetsT, weight);
      hmmS2_.bjetsM_n->Fill(nBjetsM, weight);
      hmmS2_.bjetsL_n->Fill(nBjetsL, weight);

      hmmS2_.jet1_btag->Fill(bTags[0], weight);
      hmmS2_.jet2_btag->Fill(bTags[1], weight);
      hmmS2_.jet3_btag->Fill(bTags[2], weight);
      hmmS2_.jet4_btag->Fill(bTags[3], weight);

      hmmS2_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);
    }
  }
  else if ( channel == 1 ) { // CH_ELEL
    heeS0_.bjetsT_n->Fill(nBjetsT, weight);
    heeS0_.bjetsM_n->Fill(nBjetsM, weight);
    heeS0_.bjetsL_n->Fill(nBjetsL, weight);

    heeS0_.jet1_btag->Fill(bTags[0], weight);
    heeS0_.jet2_btag->Fill(bTags[1], weight);
    heeS0_.jet3_btag->Fill(bTags[2], weight);
    heeS0_.jet4_btag->Fill(bTags[3], weight);

    heeS0_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);

    if ( nBjetsM >= 2 ) {
      heeS1_.bjetsT_n->Fill(nBjetsT, weight);
      heeS1_.bjetsM_n->Fill(nBjetsM, weight);
      heeS1_.bjetsL_n->Fill(nBjetsL, weight);

      heeS1_.jet1_btag->Fill(bTags[0], weight);
      heeS1_.jet2_btag->Fill(bTags[1], weight);
      heeS1_.jet3_btag->Fill(bTags[2], weight);
      heeS1_.jet4_btag->Fill(bTags[3], weight);

      heeS1_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);
    }

    if ( nBjetsT >= 2 ) {
      heeS2_.bjetsT_n->Fill(nBjetsT, weight);
      heeS2_.bjetsM_n->Fill(nBjetsM, weight);
      heeS2_.bjetsL_n->Fill(nBjetsL, weight);

      heeS2_.jet1_btag->Fill(bTags[0], weight);
      heeS2_.jet2_btag->Fill(bTags[1], weight);
      heeS2_.jet3_btag->Fill(bTags[2], weight);
      heeS2_.jet4_btag->Fill(bTags[3], weight);

      heeS2_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);
    }
  }
  else if ( channel == 2 ) { // CH_MUEL
    hemS0_.bjetsT_n->Fill(nBjetsT, weight);
    hemS0_.bjetsM_n->Fill(nBjetsM, weight);
    hemS0_.bjetsL_n->Fill(nBjetsL, weight);

    hemS0_.jet1_btag->Fill(bTags[0], weight);
    hemS0_.jet2_btag->Fill(bTags[1], weight);
    hemS0_.jet3_btag->Fill(bTags[2], weight);
    hemS0_.jet4_btag->Fill(bTags[3], weight);

    hemS0_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);

    if ( nBjetsM >= 2 ) {
      hemS1_.bjetsT_n->Fill(nBjetsT, weight);
      hemS1_.bjetsM_n->Fill(nBjetsM, weight);
      hemS1_.bjetsL_n->Fill(nBjetsL, weight);

      hemS1_.jet1_btag->Fill(bTags[0], weight);
      hemS1_.jet2_btag->Fill(bTags[1], weight);
      hemS1_.jet3_btag->Fill(bTags[2], weight);
      hemS1_.jet4_btag->Fill(bTags[3], weight);

      hemS1_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);
    }

    if ( nBjetsT >= 2 ) {
      hemS2_.bjetsT_n->Fill(nBjetsT, weight);
      hemS2_.bjetsM_n->Fill(nBjetsM, weight);
      hemS2_.bjetsL_n->Fill(nBjetsL, weight);

      hemS2_.jet1_btag->Fill(bTags[0], weight);
      hemS2_.jet2_btag->Fill(bTags[1], weight);
      hemS2_.jet3_btag->Fill(bTags[2], weight);
      hemS2_.jet4_btag->Fill(bTags[3], weight);

      hemS2_.jet3_btag__jet4_btag->Fill(bTags[2], bTags[3], weight);
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTBBLLAnalyzer);

