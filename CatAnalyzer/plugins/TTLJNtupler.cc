#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"

#include "CATTools/DataFormats/interface/Lepton.h"
#include "CATTools/DataFormats/interface/Jet.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "CATTools/DataFormats/interface/GenTop.h"

#include "TTree.h"
#include "TString.h"

#include <memory>

class TTLJNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  TTLJNtupler(const edm::ParameterSet& pset);
  ~TTLJNtupler() {};

  void analyze(const edm::Event& event, const edm::EventSetup&) override;

private:
  bool isMC_, isTTbar_;

  typedef std::vector<float> vfloat;
  typedef std::vector<double> vdouble;
  edm::EDGetTokenT<float> topPtWeightToken_;
  edm::EDGetTokenT<float> puWeightToken_, puUpWeightToken_, puDnWeightToken_;
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<vfloat> scaleUpWeightsToken_, scaleDnWeightsToken_;
  edm::EDGetTokenT<vfloat> pdfWeightsToken_;
  edm::EDGetTokenT<float> csvWeightToken_;
  edm::EDGetTokenT<vfloat> csvWeightSystToken_;
  edm::EDGetTokenT<int> nVertexToken_;

  edm::EDGetTokenT<cat::LeptonCollection> leptonToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<float> metPtToken_, metPhiToken_;

private:
  void clear();

  TTree* tree_;
  unsigned int b_run;
  unsigned short b_lumi;
  unsigned long long b_event;

  float b_wgt_gen;
  float b_wgt_pu, b_wgt_puUp, b_wgt_puDn;
  float b_wgt_csv;
  float b_wgt_topPt;

  std::unique_ptr<std::vector<float>> b_wgts_pdf;
  std::unique_ptr<std::vector<float>> b_wgts_scaleUp, b_wgts_scaleDn;
  std::unique_ptr<std::vector<float>> b_wgts_csvUp, b_wgts_csvDn;

  float b_event_st;
  unsigned char b_vertex_n;

  const unsigned static char kMaxNLeptons = 10;
  unsigned char b_leptons_n;
  float b_leptons_pt[kMaxNLeptons], b_leptons_eta[kMaxNLeptons], b_leptons_phi[kMaxNLeptons];
  unsigned char b_leptons_pid[kMaxNLeptons];

  float b_met_pt, b_met_phi;

  const unsigned static char kMaxNJets = 100;
  unsigned char b_jets_n;
  float b_jets_ht;
  float b_jets_pt[kMaxNJets], b_jets_eta[kMaxNJets], b_jets_phi[kMaxNJets];
  float b_jets_m[kMaxNJets];
  float b_jets_csv[kMaxNJets], b_jets_CvsL[kMaxNJets], b_jets_CvsB[kMaxNJets];

  unsigned char b_bjetsL_n, b_bjetsM_n, b_bjetsT_n;
};

using namespace std;

TTLJNtupler::TTLJNtupler(const edm::ParameterSet& pset)
{
  isMC_ = pset.getParameter<bool>("isMC");
  isTTbar_ = pset.getParameter<bool>("isTTbar");

  auto srcTag = pset.getParameter<edm::InputTag>("src");
  leptonToken_ = consumes<cat::LeptonCollection>(edm::InputTag(srcTag.label(), "leptons"));
  jetToken_ = consumes<cat::JetCollection>(edm::InputTag(srcTag.label(), "jets"));
  metPtToken_  = consumes<float>(edm::InputTag(srcTag.label(), "met"));
  metPhiToken_ = consumes<float>(edm::InputTag(srcTag.label(), "metphi"));

  nVertexToken_ = consumes<int>(pset.getParameter<edm::InputTag>("nVertex"));
  if ( isMC_ ) {
    auto puWeightLabel = pset.getParameter<edm::InputTag>("puWeight");
    puWeightToken_ = consumes<float>(puWeightLabel);
    puUpWeightToken_ = consumes<float>(edm::InputTag(puWeightLabel.label(), "up"));
    puDnWeightToken_ = consumes<float>(edm::InputTag(puWeightLabel.label(), "dn"));

    auto genWeightLabel = pset.getParameter<edm::InputTag>("genWeight");
    genWeightToken_ = consumes<float>(genWeightLabel);
    scaleUpWeightsToken_ = consumes<vfloat>(edm::InputTag(genWeightLabel.label(), "scaleup"));
    scaleDnWeightsToken_ = consumes<vfloat>(edm::InputTag(genWeightLabel.label(), "scaledown"));
    pdfWeightsToken_ = consumes<vfloat>(edm::InputTag(genWeightLabel.label(), "pdf"));

    auto csvWeightLabel = pset.getParameter<edm::InputTag>("csvWeight");
    csvWeightToken_ = consumes<float>(edm::InputTag(csvWeightLabel.label()));
    csvWeightSystToken_ = consumes<vfloat>(edm::InputTag(csvWeightLabel.label(), "syst"));

    topPtWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("topPtWeight"));
  }

  b_wgts_scaleUp.reset(new std::vector<float>());
  b_wgts_scaleDn.reset(new std::vector<float>());
  b_wgts_pdf.reset(new std::vector<float>());
  b_wgts_csvUp.reset(new std::vector<float>());
  b_wgts_csvDn.reset(new std::vector<float>());

  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("events", "TTLJ event tree");

  tree_->Branch("run", &b_run, "run/i"); // run number in 10^7
  tree_->Branch("lumi", &b_lumi, "lumi/s"); // enough with 65535
  tree_->Branch("event", &b_event, "event/l"); // event number with 64bit int

  if ( isMC_ ) {
    // event weights
    tree_->Branch("wgt_gen", &b_wgt_gen, "wgt_gen/F"); // central gen weight
    tree_->Branch("wgt_pu", &b_wgt_pu, "wgt_pu/F"); // pileup weight
    tree_->Branch("wgt_csv", &b_wgt_csv, "wgt_csv/F"); // CSV weight
    if ( isTTbar_ ) {
      tree_->Branch("wgt_topPt", &b_wgt_topPt, "wgt_topPt/F");// top pt weight
    }

    tree_->Branch("wgt_puUp", &b_wgt_puUp, "wgt_puUp/F"); // pileup weight
    tree_->Branch("wgt_puDn", &b_wgt_puDn, "wgt_puDn/F"); // pileup weight

    tree_->Branch("wgts_scaleUp", "std::vector<float>", &*b_wgts_scaleUp);
    tree_->Branch("wgts_scaleDn", "std::vector<float>", &*b_wgts_scaleDn);
    tree_->Branch("wgts_pdf", "std::vector<float>", &*b_wgts_pdf);

    tree_->Branch("wgts_csvUp", "std::vector<float>", &*b_wgts_csvUp); // CSV weight
    tree_->Branch("wgts_csvDn", "std::vector<float>", &*b_wgts_csvDn); // CSV weight
  }

  // Event-wise observables
  tree_->Branch("event_st", &b_event_st, "event_st/F");
  tree_->Branch("vertex_n", &b_vertex_n, "vertex_n/b"); // enough with 255

  // Leptons
  tree_->Branch("leptons_n", &b_leptons_n, "leptons_n/b");
  tree_->Branch("leptons_pt" , &b_leptons_pt , "leptons_pt[leptons_n]/F");
  tree_->Branch("leptons_eta", &b_leptons_eta, "leptons_eta[leptons_n]/F");
  tree_->Branch("leptons_phi", &b_leptons_phi, "leptons_phi[leptons_n]/F");
  tree_->Branch("leptons_pid", &b_leptons_pid, "leptons_pid[leptons_n]/B"); // +-11 or +-13. enough with += 127

  // MET
  tree_->Branch("met_pt" , &b_met_pt , "met_pt/F");
  tree_->Branch("met_phi", &b_met_phi, "met_phi/F");

  // Jets
  tree_->Branch("jets_n", &b_jets_n, "jets_n/b"); // enough with 255
  tree_->Branch("jets_ht", &b_jets_ht, "jets_ht/F");
  tree_->Branch("jets_pt" , &b_jets_pt , "jets_pt[jets_n]/F");
  tree_->Branch("jets_eta", &b_jets_eta, "jets_eta[jets_n]/F");
  tree_->Branch("jets_phi", &b_jets_phi, "jets_phi[jets_n]/F");
  tree_->Branch("jets_m"  , &b_jets_m  , "jets_m[jets_n]/F");
  tree_->Branch("jets_csv" , &b_jets_csv , "jets_csv[jets_n]/F");
  tree_->Branch("jets_CvsL", &b_jets_CvsL, "jets_CvsL[jets_n]/F");
  tree_->Branch("jets_CvsB", &b_jets_CvsB, "jets_CvsB[jets_n]/F");

  // bJets
  tree_->Branch("bjetsL_n", &b_bjetsL_n, "bjetsL_n/b"); // enough with 255
  tree_->Branch("bjetsM_n", &b_bjetsM_n, "bjetsM_n/b"); // enough with 255
  tree_->Branch("bjetsT_n", &b_bjetsT_n, "bjetsT_n/b"); // enough with 255

  // Composite objects
}

void TTLJNtupler::analyze(const edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<float> fHandle;
  edm::Handle<int> iHandle;
  edm::Handle<std::vector<float>> vfHandle;

  clear();
  b_run   = event.id().run();
  b_lumi  = event.id().luminosityBlock();
  b_event = event.id().event();

  if ( isMC_ ) {
    if ( isTTbar_ and event.getByToken(topPtWeightToken_, fHandle) ) {
      b_wgt_topPt = *fHandle;
    }

    event.getByToken(puWeightToken_, fHandle);
    b_wgt_pu = *fHandle;

    event.getByToken(puUpWeightToken_, fHandle);
    b_wgt_puUp = *fHandle;

    event.getByToken(puDnWeightToken_, fHandle);
    b_wgt_puDn = *fHandle;

    event.getByToken(genWeightToken_, fHandle);
    b_wgt_gen = *fHandle;

    if ( event.getByToken(scaleUpWeightsToken_, vfHandle) ) {
      for ( auto x : *vfHandle ) b_wgts_scaleUp->push_back(x);
    }
    if ( event.getByToken(scaleDnWeightsToken_, vfHandle) ) {
      for ( auto x : *vfHandle ) b_wgts_scaleDn->push_back(x);
    }

    if ( event.getByToken(pdfWeightsToken_, vfHandle) ) {
      for ( auto x : *vfHandle ) b_wgts_pdf->push_back(x);
    }

    event.getByToken(csvWeightToken_, fHandle);
    b_wgt_csv = *fHandle;
    event.getByToken(csvWeightSystToken_, vfHandle);
    for ( int i=0, n=vfHandle->size(); i<n; i+=2 ) {
      b_wgts_csvUp->push_back(vfHandle->at(i+0));
      b_wgts_csvDn->push_back(vfHandle->at(i+1));
    }
  }

  event.getByToken(nVertexToken_, iHandle);
  b_vertex_n = *iHandle;

  event.getByToken(metPtToken_, fHandle);
  b_met_pt  = *fHandle;
  b_event_st += b_met_pt;

  event.getByToken(metPhiToken_, fHandle);
  b_met_phi = *fHandle;

  edm::Handle<cat::LeptonCollection> leptonHandle;
  event.getByToken(leptonToken_, leptonHandle);
  if ( !leptonHandle->empty() ) {
    b_leptons_n = 1;
    const auto& lepton1 = leptonHandle->at(0);
    b_leptons_pt[0] = lepton1.pt();
    b_leptons_eta[0] = lepton1.eta();
    b_leptons_phi[0] = lepton1.phi();
    b_leptons_pid[0] = lepton1.pdgId();

    b_event_st += lepton1.pt();
  }

  edm::Handle<cat::JetCollection> jetHandle;
  event.getByToken(jetToken_, jetHandle);
  b_jets_n = jetHandle->size();
  for ( int i=0, n=jetHandle->size(); i<n; ++i ) {
    const auto& jet = jetHandle->at(i);
    if ( i < kMaxNJets ) {
      b_jets_pt[i] = jet.pt();
      b_jets_eta[i] = jet.eta();
      b_jets_phi[i] = jet.phi();
      b_jets_m[i] = jet.mass();
      b_jets_csv[i] = jet.bDiscriminator(cat::BTAG_CSVv2);
      b_jets_CvsL[i] = jet.bDiscriminator(cat::CTAG_CvsL);
      b_jets_CvsB[i] = jet.bDiscriminator(cat::CTAG_CvsB);
    }
    b_jets_ht += jet.pt();

    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2L ) ++b_bjetsL_n;
    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2M ) ++b_bjetsM_n;
    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2T ) ++b_bjetsT_n;
  }
  b_event_st += b_jets_ht;

  tree_->Fill();
}

void TTLJNtupler::clear()
{
  b_run = b_lumi = b_event = 0;

  b_wgt_gen = 1.;
  b_wgt_pu = b_wgt_puUp = b_wgt_puDn = 1.;
  b_wgt_csv = 1.;
  b_wgt_topPt = 1.;

  b_wgts_pdf->clear();
  b_wgts_scaleUp->clear();
  b_wgts_scaleDn->clear();
  b_wgts_csvUp->clear();
  b_wgts_csvDn->clear();

  b_event_st = -10;
  b_vertex_n = 0;
  b_leptons_n = 0;
  b_leptons_pt[0] = b_leptons_eta[0] = b_leptons_phi[0] = -10;
  b_leptons_pid[0] = 0;
  b_met_pt = b_met_phi = -10;
  b_jets_n = 0;
  b_jets_ht = -10;
  for ( int i=0; i<kMaxNJets; ++i ) {
    b_jets_pt[i] = b_jets_eta[i] = b_jets_phi[i] = -10;
    b_jets_m[i] = -10;
    b_jets_csv[i] = b_jets_CvsL[i] = b_jets_CvsB[i] = -10;
  }
  b_bjetsL_n = b_bjetsM_n = b_bjetsT_n = 0;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLJNtupler);
