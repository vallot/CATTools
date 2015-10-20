#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/PtEtaPhiM4D.h"

#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

class ColorCoherenceAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ColorCoherenceAnalyzer(const edm::ParameterSet&);
  ~ColorCoherenceAnalyzer() {};

  enum sys_e {sys_nom, sys_jes_u, sys_jes_d, sys_jer_u, sys_jer_d, sys_jar, nsys_e};

private:

  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void resetBr();
  vector<TLorentzVector> selectJets(const cat::JetCollection& jets, sys_e sys) const;
  TLorentzVector sysJet(const cat::Jet&, sys_e sys) const;
  TLorentzVector jarJet(const cat::Jet& jet) const;

  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<int>   vtxToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  vector<TTree*> ttree_;

  int b_nVtx, b_nJet, b_hlt_40_pass, b_hlt_60_pass, b_hlt_80_pass, b_hlt_140_pass, b_hlt_320_pass, b_hlt_400_pass, b_hlt_450_pass, b_hlt_500_pass;
  float b_beta, b_del_eta, b_del_phi, b_del_r;
  float b_del_r12, b_raw_mass;
  float b_jet1_pt, b_jet1_eta, b_jet1_phi;
  float b_jet2_pt, b_jet2_eta, b_jet2_phi;
  float b_jet3_pt, b_jet3_eta, b_jet3_phi;
  float b_met, b_metSig;

  bool runOnMC_;
  float b_gjet1_pt, b_gjet1_eta, b_gjet1_phi;
  float b_gjet2_pt, b_gjet2_eta, b_gjet2_phi;
  float b_gjet3_pt, b_gjet3_eta, b_gjet3_phi;
  float b_gbeta, b_gdel_eta, b_gdel_phi, b_gdel_r;

};

ColorCoherenceAnalyzer::ColorCoherenceAnalyzer(const edm::ParameterSet& iConfig)
{
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<int>(iConfig.getParameter<edm::InputTag>("vtx"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  const std::vector<std::string> sys_name = {"nom", "jes_u", "jes_d", "jer_u", "jer_d", "jar"};
  for (auto sys : sys_name) {
    ttree_.push_back(fs->make<TTree>(sys.c_str(), (string("color cohernece systematic errors : ")+sys).c_str()));
    auto tr = ttree_.back();

    tr->Branch("nvtx", &b_nVtx, "nvtx/I");
    tr->Branch("njet", &b_nJet, "njet/I");
    tr->Branch("hlt_40_pass", &b_hlt_40_pass, "hlt_40_pass/I");
    tr->Branch("hlt_60_pass", &b_hlt_60_pass, "hlt_60_pass/I");
    tr->Branch("hlt_80_pass", &b_hlt_80_pass, "hlt_80_pass/I");

    tr->Branch("hlt_140_pass", &b_hlt_140_pass, "hlt_140_pass/I");
    tr->Branch("hlt_320_pass", &b_hlt_320_pass, "hlt_320_pass/I");
    tr->Branch("hlt_400_pass", &b_hlt_400_pass, "hlt_400_pass/I");
    tr->Branch("hlt_450_pass", &b_hlt_450_pass, "hlt_450_pass/I");
    tr->Branch("hlt_500_pass", &b_hlt_500_pass, "hlt_500_pass/I");
    tr->Branch("beta", &b_beta, "beta/F");
    tr->Branch("del_eta", &b_del_eta, "del_eta/F");
    tr->Branch("del_phi", &b_del_phi, "del_phi/F");
    tr->Branch("del_r", &b_del_r, "del_r/F");
    tr->Branch("del_r12", &b_del_r12, "del_r12/F");
    tr->Branch("raw_mass", &b_raw_mass, "raw_mass/F");
    tr->Branch("jet1_pt", &b_jet1_pt, "jet1_pt/F");
    tr->Branch("jet1_eta", &b_jet1_eta, "jet1_eta/F");
    tr->Branch("jet1_phi", &b_jet1_phi, "jet1_phi/F");
    tr->Branch("jet2_pt", &b_jet2_pt, "jet2_pt/F");
    tr->Branch("jet2_eta", &b_jet2_eta, "jet2_eta/F");
    tr->Branch("jet2_phi", &b_jet2_phi, "jet2_phi/F");
    tr->Branch("jet3_pt", &b_jet3_pt, "jet3_pt/F");
    tr->Branch("jet3_eta", &b_jet3_eta, "jet3_eta/F");
    tr->Branch("jet3_phi", &b_jet3_phi, "jet3_phi/F");
    tr->Branch("met", &b_met, "met/F");
    tr->Branch("metSig", &b_metSig, "metSig/F");
    tr->Branch("gjet1_pt", &b_gjet1_pt, "gjet1_pt/F");
    tr->Branch("gjet1_eta", &b_gjet1_eta, "gjet1_eta/F");
    tr->Branch("gjet1_phi", &b_gjet1_phi, "gjet1_phi/F");
    tr->Branch("gjet2_pt", &b_gjet2_pt, "gjet2_pt/F");
    tr->Branch("gjet2_eta", &b_gjet2_eta, "gjet2_eta/F");
    tr->Branch("gjet2_phi", &b_gjet2_phi, "gjet2_phi/F");
    tr->Branch("gjet3_pt", &b_gjet3_pt, "gjet3_pt/F");
    tr->Branch("gjet3_eta", &b_gjet3_eta, "gjet3_eta/F");
    tr->Branch("gjet3_phi", &b_gjet3_phi, "gjet3_phi/F");
    tr->Branch("gbeta", &b_gbeta, "gbeta/F");
    tr->Branch("gdel_eta", &b_gdel_eta, "gdel_eta/F");
    tr->Branch("gdel_phi", &b_gdel_phi, "gdel_phi/F");
    tr->Branch("gdel_r", &b_gdel_r, "gdel_r/F");
  }

}

void ColorCoherenceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  edm::Handle<cat::JetCollection> jets;
  if (!iEvent.getByToken(jetToken_, jets)) return;

  edm::Handle<cat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);

  edm::Handle<int> vtx;
  iEvent.getByToken(vtxToken_, vtx);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  cat::AnalysisHelper trigHelper = cat::AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  for (int sys = 0; sys < nsys_e; ++sys){
    resetBr();

    vector<TLorentzVector>&& seljets = selectJets(*(jets.product()), (sys_e)sys);
    if (seljets.size() < 3) return;

    sort(seljets.begin(), seljets.end(), cat::GtByTLVPt);

    b_nVtx = *vtx;
    b_nJet = seljets.size();

    int hlt_count = 0;
    //trigHelper.listFiredTriggers();
    if(trigHelper.triggerFired("HLT_PAJet40_NoJetID_v")) {b_hlt_40_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PAJet60_NoJetID_v")) {b_hlt_60_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PAJet80_NoJetID_v")) {b_hlt_80_pass = 1; hlt_count++;}

    if(trigHelper.triggerFired("HLT_PFJet80_v")) {b_hlt_80_pass = 1; hlt_count++;}// dont u use pfjet80 too?
    if(trigHelper.triggerFired("HLT_PFJet140_v")) {b_hlt_140_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet320_v")) {b_hlt_320_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet400_v")) {b_hlt_400_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet450_v")) {b_hlt_450_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet500_v")) {b_hlt_500_pass = 1; hlt_count++;}
    //if (hlt_count < 1) return;

    b_del_eta = copysign(1.0, seljets[1].Eta())*(seljets[2].Eta() - seljets[1].Eta());
    b_del_phi = reco::deltaPhi(seljets[2].Phi(), seljets[1].Phi());
    b_beta = atan2(b_del_phi, b_del_eta);
    b_del_r = reco::deltaR(seljets[2].Eta(), seljets[2].Phi(), seljets[1].Eta(), seljets[1].Phi());
    b_del_r12 = reco::deltaR(seljets[1].Eta(), seljets[1].Phi(), seljets[0].Eta(), seljets[0].Phi());
    b_raw_mass = (seljets[0] + seljets[1]).M();

    b_jet1_pt = seljets[0].Pt(); b_jet1_eta = seljets[0].Eta(); b_jet1_phi = seljets[0].Phi();
    b_jet2_pt = seljets[1].Pt(); b_jet2_eta = seljets[1].Eta(); b_jet2_phi = seljets[1].Phi();
    b_jet3_pt = seljets[2].Pt(); b_jet3_eta = seljets[2].Eta(); b_jet3_phi = seljets[2].Phi();

    // TODO: MET uncertainty due to jet energy scale to be added!!!
    b_met = mets->begin()->et();
    b_metSig = mets->begin()->et()/mets->begin()->sumEt();
    ttree_[sys]->Fill();

  }
}

vector<TLorentzVector> ColorCoherenceAnalyzer::selectJets(const cat::JetCollection& jets, sys_e sys) const
{
  vector<TLorentzVector> seljets;
  for (auto jet : jets) {
    if (!jet.LooseId()) continue;
    if (jet.pileupJetId() <0.9) continue;
    const TLorentzVector&& newjet = sysJet(jet, sys);
    if (newjet.Pt() <= 20.) continue;

    seljets.push_back(newjet);
  }
  return seljets;
}

TLorentzVector ColorCoherenceAnalyzer::sysJet(const cat::Jet& jet, sys_e sys) const
{
  if (sys == sys_nom) return jet.tlv();
  if (sys == sys_jes_u) return jet.tlv()*jet.shiftedEnUp();
  if (sys == sys_jes_d) return jet.tlv()*jet.shiftedEnDown();
  if (sys == sys_jer_u) return jet.tlv()*jet.smearedResUp();
  if (sys == sys_jer_d) return jet.tlv()*jet.smearedResDown();
  if (sys == sys_jar) return jarJet(jet);

  return jet.tlv();
}

TLorentzVector ColorCoherenceAnalyzer::jarJet(const cat::Jet& jet) const
{
  // temp... currently doing NOTHING
  TLorentzVector sJet;
  sJet.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.p4().M());
  return sJet;
}

void ColorCoherenceAnalyzer::resetBr()
{
  b_nVtx = -99; b_nJet = -99; b_hlt_40_pass = -99; b_hlt_60_pass = -99; b_hlt_80_pass = -99; b_hlt_140_pass = -99;  b_hlt_320_pass = -99; b_hlt_400_pass = -99; b_hlt_450_pass = -99; b_hlt_500_pass = -99;
  b_beta = -99; b_del_eta = -99; b_del_phi = -99; b_del_r = -99;
  b_del_r12 = -99; b_raw_mass = -99;
  b_jet1_pt = -99; b_jet1_eta = -99; b_jet1_phi = -99;
  b_jet2_pt = -99; b_jet2_eta = -99; b_jet2_phi = -99;
  b_jet3_pt = -99; b_jet3_eta = -99; b_jet3_phi = -99;
  b_met = -99;

  b_gbeta = -99; b_gdel_eta = -99; b_gdel_phi = -99; b_gdel_r = -99;
  b_gjet1_pt = -99; b_gjet1_eta = -99; b_gjet1_phi = -99;
  b_gjet2_pt = -99; b_gjet2_eta = -99; b_gjet2_phi = -99;
  b_gjet3_pt = -99; b_gjet3_eta = -99; b_gjet3_phi = -99;
}

DEFINE_FWK_MODULE(ColorCoherenceAnalyzer);
