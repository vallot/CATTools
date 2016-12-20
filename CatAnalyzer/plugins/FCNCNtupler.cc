#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"

#include "CATTools/DataFormats/interface/Lepton.h"
#include "CATTools/DataFormats/interface/Jet.h"
/*
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenWeights.h"

#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"

#include "TH1.h"
*/

#include "TTree.h"
#include "TString.h"

#include <memory>

class FCNCNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  FCNCNtupler(const edm::ParameterSet& pset);
  ~FCNCNtupler() {};

  void analyze(const edm::Event& event, const edm::EventSetup&) override;

private:
  bool isMC_; // set to true by default. just let the framework to decide behaviour

  typedef std::vector<float> vfloat;
  typedef std::vector<double> vdouble;
  edm::EDGetTokenT<float> puWeightToken_, genWeightToken_;
  edm::EDGetTokenT<vfloat> scaleWeightsToken_, pdfWeightsToken_;
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

  float b_weight_gen, b_weight_pu;
  float b_weight_csv;

  float b_event_st;
  unsigned char b_vertex_n;

  float b_lepton1_pt, b_lepton1_eta, b_lepton1_phi;
  unsigned char b_lepton1_pid;

  float b_met_pt, b_met_phi;

  const unsigned static char kMaxNJets = 6;
  unsigned char b_jets_n;
  float b_jets_ht;
  float b_jets_pt[kMaxNJets], b_jets_eta[kMaxNJets], b_jets_phi[kMaxNJets];
  float b_jets_m[kMaxNJets];
  float b_jets_btag[kMaxNJets], b_jets_ctagL[kMaxNJets], b_jets_ctagB[kMaxNJets];

  unsigned char b_bjets_n;
};

using namespace std;

FCNCNtupler::FCNCNtupler(const edm::ParameterSet& pset)
{
  isMC_ = pset.getParameter<bool>("isMC");

  auto srcTag = pset.getParameter<edm::InputTag>("src");
  leptonToken_ = consumes<cat::LeptonCollection>(edm::InputTag(srcTag.label(), "leptons"));
  jetToken_ = consumes<cat::JetCollection>(edm::InputTag(srcTag.label(), "jets"));
  metPtToken_  = consumes<float>(edm::InputTag(srcTag.label(), "met"));
  metPhiToken_ = consumes<float>(edm::InputTag(srcTag.label(), "metphi"));

  nVertexToken_ = consumes<int>(pset.getParameter<edm::InputTag>("nVertex"));
  if ( isMC_ ) {
    puWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("puWeight"));
    genWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("genWeight"));
    scaleWeightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("scaleWeights"));
    pdfWeightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("pdfWeights"));
    csvWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("csvWeight"));
    csvWeightSystToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("csvWeightSyst"));
  }

  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("events", "FCNC event tree");

  tree_->Branch("run", &b_run, "run/i"); // run number in 10^7
  tree_->Branch("lumi", &b_lumi, "lumi/s"); // enough with 65535
  tree_->Branch("event", &b_event, "event/l"); // event number with 64bit int

  // event weights
  tree_->Branch("weight_gen", &b_weight_gen, "weight_gen/F"); // central gen weight
  tree_->Branch("weight_pu", &b_weight_pu, "weight_pu/F"); // pileup weight
  tree_->Branch("weight_csv", &b_weight_csv, "weight_csv/F"); // CSV weight

  // Event-wise observables
  tree_->Branch("event_st", &b_event_st, "event_st/F");
  tree_->Branch("vertex_n", &b_vertex_n, "vertex_n/b"); // enough with 255

  // Leptons
  tree_->Branch("lepton1_pt" , &b_lepton1_pt , "lepton1_pt/F");
  tree_->Branch("lepton1_eta", &b_lepton1_eta, "lepton1_eta/F");
  tree_->Branch("lepton1_phi", &b_lepton1_phi, "lepton1_phi/F");
  tree_->Branch("lepton1_pid", &b_lepton1_pid, "lepton1_pid/B"); // +-11 or +-13. enough with += 127

  // MET
  tree_->Branch("met_pt" , &b_met_pt , "met_pt/F");
  tree_->Branch("met_phi", &b_met_phi, "met_phi/F");

  // Jets
  tree_->Branch("jets_n", &b_jets_n, "jets_n/b"); // enough with 255
  tree_->Branch("jets_ht", &b_jets_ht, "jets_ht/F");
  for ( int i=0; i<kMaxNJets; ++i ) {
    const string name = Form("jet%d_", i);
    tree_->Branch((name+"_pt" ).c_str(), &b_jets_pt[i] , (name+"_pt/F" ).c_str());
    tree_->Branch((name+"_eta").c_str(), &b_jets_eta[i], (name+"_eta/F").c_str());
    tree_->Branch((name+"_phi").c_str(), &b_jets_phi[i], (name+"_phi/F").c_str());
    tree_->Branch((name+"_m"  ).c_str(), &b_jets_m[i]  , (name+"_m/F"  ).c_str());
    tree_->Branch((name+"_btag").c_str(), &b_jets_btag[i], (name+"_btag/F").c_str());
    tree_->Branch((name+"_ctagL").c_str(), &b_jets_ctagL[i], (name+"_ctagL/F").c_str());
    tree_->Branch((name+"_ctagB").c_str(), &b_jets_ctagB[i], (name+"_ctagB/F").c_str());
  }

  // bJets
  tree_->Branch("bjets_n", &b_bjets_n, "bjets_n/b"); // enough with 255

  // Composite objects
}

void FCNCNtupler::analyze(const edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<float> fHandle;
  edm::Handle<int> iHandle;

  clear();
  b_run   = event.id().run();
  b_lumi  = event.id().luminosityBlock();
  b_event = event.id().event();

  if ( isMC_ ) {
    event.getByToken(puWeightToken_, fHandle);
    b_weight_pu = *fHandle;

    event.getByToken(genWeightToken_, fHandle);
    b_weight_gen = *fHandle;

    event.getByToken(csvWeightToken_, fHandle);
    b_weight_gen = *fHandle;
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
    const auto& lepton1 = leptonHandle->at(0);
    b_lepton1_pt = lepton1.pt();
    b_lepton1_eta = lepton1.eta();
    b_lepton1_phi = lepton1.phi();
    b_lepton1_pid = lepton1.pdgId();

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
      b_jets_btag[i] = jet.bDiscriminator(cat::BTAG_CSVv2);
      b_jets_ctagL[i] = jet.bDiscriminator(cat::CTAG_CvsL);
      b_jets_ctagB[i] = jet.bDiscriminator(cat::CTAG_CvsB);
    }
    b_jets_ht += jet.pt();

    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2M ) ++b_bjets_n;
  }
  b_event_st += b_jets_ht;

  tree_->Fill();
}

void FCNCNtupler::clear()
{
  b_run = b_lumi = b_event = 0;
  b_weight_gen = b_weight_pu = b_weight_csv = 1.;
  b_event_st = -1e4;
  b_vertex_n = 0;
  b_lepton1_pt = b_lepton1_eta = b_lepton1_phi = -1e4;
  b_lepton1_pid = 0;
  b_met_pt = b_met_phi = -1e4;
  b_jets_n = 0;
  b_jets_ht = -1e4;
  for ( int i=0; i<kMaxNJets; ++i ) {
    b_jets_pt[i] = b_jets_eta[i] = b_jets_phi[i] = -1;
    b_jets_m[i] = -1;
    b_jets_btag[i] = b_jets_ctagL[i] = b_jets_ctagB[i] = -1e4;
  }
  b_bjets_n = 0;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FCNCNtupler);
