#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

using namespace std;
using namespace cat;

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

  typedef TH1D* H1;
  H1 hee_bjetsT_n_, hee_bjetsM_n_, hee_bjetsL_n_;
  H1 hee_bsortedjet1_btag_, hee_bsortedjet2_btag_, hee_bsortedjet3_btag_, hee_bsortedjet4_btag_;

  H1 hmm_bjetsT_n_, hmm_bjetsM_n_, hmm_bjetsL_n_;
  H1 hmm_bsortedjet1_btag_, hmm_bsortedjet2_btag_, hmm_bsortedjet3_btag_, hmm_bsortedjet4_btag_;

  H1 hem_bjetsT_n_, hem_bjetsM_n_, hem_bjetsL_n_;
  H1 hem_bsortedjet1_btag_, hem_bsortedjet2_btag_, hem_bsortedjet3_btag_, hem_bsortedjet4_btag_;
};

TTBBLLAnalyzer::TTBBLLAnalyzer(const edm::ParameterSet& pset)
{
  channelToken_ = consumes<int>(edm::InputTag("ttll", "channel"));
  weightToken_ = consumes<float>(edm::InputTag("ttll", "weight"));
  leptonsToken_ = consumes<cat::LeptonCollection>(edm::InputTag("ttll", "leptons"));
  jetsToken_ = consumes<cat::JetCollection>(edm::InputTag("ttll", "jets"));
  metToken_ = consumes<float>(edm::InputTag("ttll", "met"));
  metphiToken_ = consumes<float>(edm::InputTag("ttll", "metphi"));

  edm::Service<TFileService> fs;
  auto diree = fs->mkdir("ee");
  hee_bjetsT_n_ = diree.make<TH1D>("bjetsT_n", "Tight b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hee_bjetsM_n_ = diree.make<TH1D>("bjetsM_n", "Medium b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hee_bjetsL_n_ = diree.make<TH1D>("bjetsL_n", "Loose b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hee_bsortedjet1_btag_ = diree.make<TH1D>("bsortedJet1_btag", "1st b discriminator;B discriminator", 100, 0, 1);
  hee_bsortedjet2_btag_ = diree.make<TH1D>("bsortedJet1_btag", "2nd b discriminator;B discriminator", 100, 0, 1);
  hee_bsortedjet3_btag_ = diree.make<TH1D>("bsortedJet1_btag", "3rd b discriminator;B discriminator", 100, 0, 1);
  hee_bsortedjet4_btag_ = diree.make<TH1D>("bsortedJet1_btag", "4th b discriminator;B discriminator", 100, 0, 1);

  auto dirmm = fs->mkdir("mm");
  hmm_bjetsT_n_ = dirmm.make<TH1D>("bjetsT_n", "Tight b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hmm_bjetsM_n_ = dirmm.make<TH1D>("bjetsM_n", "Medium b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hmm_bjetsL_n_ = dirmm.make<TH1D>("bjetsL_n", "Loose b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hmm_bsortedjet1_btag_ = dirmm.make<TH1D>("bsortedJet1_btag", "1st b discriminator;B discriminator", 100, 0, 1);
  hmm_bsortedjet2_btag_ = dirmm.make<TH1D>("bsortedJet1_btag", "2nd b discriminator;B discriminator", 100, 0, 1);
  hmm_bsortedjet3_btag_ = dirmm.make<TH1D>("bsortedJet1_btag", "3rd b discriminator;B discriminator", 100, 0, 1);
  hmm_bsortedjet4_btag_ = dirmm.make<TH1D>("bsortedJet1_btag", "4th b discriminator;B discriminator", 100, 0, 1);

  auto direm = fs->mkdir("em");
  hem_bjetsT_n_ = direm.make<TH1D>("bjetsT_n", "Tight b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hem_bjetsM_n_ = direm.make<TH1D>("bjetsM_n", "Medium b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hem_bjetsL_n_ = direm.make<TH1D>("bjetsL_n", "Loose b jet multiplicity;B jet multiplicity", 10, 0, 10);
  hem_bsortedjet1_btag_ = direm.make<TH1D>("bsortedJet1_btag", "1st b discriminator;B discriminator", 100, 0, 1);
  hem_bsortedjet2_btag_ = direm.make<TH1D>("bsortedJet1_btag", "2nd b discriminator;B discriminator", 100, 0, 1);
  hem_bsortedjet3_btag_ = direm.make<TH1D>("bsortedJet1_btag", "3rd b discriminator;B discriminator", 100, 0, 1);
  hem_bsortedjet4_btag_ = direm.make<TH1D>("bsortedJet1_btag", "4th b discriminator;B discriminator", 100, 0, 1);
  
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
    hmm_bjetsT_n_->Fill(nBjetsT, weight);
    hmm_bjetsM_n_->Fill(nBjetsM, weight);
    hmm_bjetsL_n_->Fill(nBjetsL, weight);

    hmm_bsortedjet1_btag_->Fill(bTags[0], weight);
    hmm_bsortedjet2_btag_->Fill(bTags[1], weight);
    hmm_bsortedjet3_btag_->Fill(bTags[2], weight);
    hmm_bsortedjet4_btag_->Fill(bTags[3], weight);
  }
  else if ( channel == 1 ) { // CH_ELEL
    hee_bjetsT_n_->Fill(nBjetsT, weight);
    hee_bjetsM_n_->Fill(nBjetsM, weight);
    hee_bjetsL_n_->Fill(nBjetsL, weight);

    hee_bsortedjet1_btag_->Fill(bTags[0], weight);
    hee_bsortedjet2_btag_->Fill(bTags[1], weight);
    hee_bsortedjet3_btag_->Fill(bTags[2], weight);
    hee_bsortedjet4_btag_->Fill(bTags[3], weight);
  }
  else if ( channel == 2 ) { // CH_MUEL
    hem_bjetsT_n_->Fill(nBjetsT, weight);
    hem_bjetsM_n_->Fill(nBjetsM, weight);
    hem_bjetsL_n_->Fill(nBjetsL, weight);

    hem_bsortedjet1_btag_->Fill(bTags[0], weight);
    hem_bsortedjet2_btag_->Fill(bTags[1], weight);
    hem_bsortedjet3_btag_->Fill(bTags[2], weight);
    hem_bsortedjet4_btag_->Fill(bTags[3], weight);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTBBLLAnalyzer);

