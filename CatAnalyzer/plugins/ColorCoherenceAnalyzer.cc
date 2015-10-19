#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "CATTools/CatAnalyzer/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/PtEtaPhiM4D.h"

#define _USE_MATH_DEFINES
#include <cmath>
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > simpleJet;
typedef GreaterByPt<simpleJet> GtByPt;
using namespace std;

class ColorCoherenceAnalyzer : public edm::EDAnalyzer{
public:
  explicit ColorCoherenceAnalyzer(const edm::ParameterSet&); 
  ~ColorCoherenceAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:

  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  vector<cat::Jet> selectJets(const edm::View<cat::Jet>* jets);
  vector<simpleJet> sortJets(vector<cat::Jet> seljets,  TString sys_name);
  simpleJet sysJet(cat::Jet, TString sys_name);
  simpleJet jarJet(cat::Jet jet);
  vector<float> calBeta(vector<simpleJet> seljets);
  void resetBr();

  edm::EDGetTokenT<edm::View<cat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> >      metToken_;
  edm::EDGetTokenT<int>   vtxToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  vector<TTree*> ttree_;
  vector<TString> sys_e;

  int b_nVtx, b_nJet, b_hlt_40_pass, b_hlt_60_pass, b_hlt_80_pass, b_hlt_140_pass, b_hlt_320_pass, b_hlt_400_pass, b_hlt_450_pass, b_hlt_500_pass;
  float b_beta, b_del_eta, b_del_phi, b_del_r;
  float b_del_r12, b_raw_mass;
  float b_jet1_pt, b_jet1_eta, b_jet1_phi;
  float b_jet2_pt, b_jet2_eta, b_jet2_phi;
  float b_jet3_pt, b_jet3_eta, b_jet3_phi;
  float b_met;

  bool runOnMC_;
  float b_gjet1_pt, b_gjet1_eta, b_gjet1_phi;
  float b_gjet2_pt, b_gjet2_eta, b_gjet2_phi;
  float b_gjet3_pt, b_gjet3_eta, b_gjet3_phi;
  float b_gbeta, b_gdel_eta, b_gdel_phi, b_gdel_r;

};

ColorCoherenceAnalyzer::ColorCoherenceAnalyzer(const edm::ParameterSet& iConfig)
{
  jetToken_  = consumes<edm::View<cat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<edm::View<cat::MET> >(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<int>(iConfig.getParameter<edm::InputTag>("vtx"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  edm::Service<TFileService> fs;
  sys_e.push_back("pt");
  sys_e.push_back("eup");
  sys_e.push_back("edown");
  sys_e.push_back("esup");
  sys_e.push_back("esdown");
  sys_e.push_back("jar");
  for (auto sys : sys_e) {
    ttree_.push_back(fs->make<TTree>(sys, "color cohernece systematic errors : "+sys));
  }
  for (auto tr : ttree_) {
    TString sys = TString(tr->GetName())+"_";

    tr->Branch(sys+"nvtx", &b_nVtx, sys+"nvtx/I");
    tr->Branch(sys+"njet", &b_nJet, sys+"njet/I");
    tr->Branch(sys+"hlt_40_pass", &b_hlt_40_pass, sys+"hlt_40_pass/I");
    tr->Branch(sys+"hlt_60_pass", &b_hlt_60_pass, sys+"hlt_60_pass/I");
    tr->Branch(sys+"hlt_80_pass", &b_hlt_80_pass, sys+"hlt_80_pass/I");
    tr->Branch(sys+"hlt_140_pass", &b_hlt_140_pass, sys+"hlt_140_pass/I");
    tr->Branch(sys+"hlt_320_pass", &b_hlt_320_pass, sys+"hlt_320_pass/I");
    tr->Branch(sys+"hlt_400_pass", &b_hlt_400_pass, sys+"hlt_400_pass/I");
    tr->Branch(sys+"hlt_450_pass", &b_hlt_450_pass, sys+"hlt_450_pass/I");
    tr->Branch(sys+"hlt_500_pass", &b_hlt_500_pass, sys+"hlt_500_pass/I");
    tr->Branch(sys+"beta", &b_beta, sys+"beta/F");
    tr->Branch(sys+"del_eta", &b_del_eta, sys+"del_eta/F");
    tr->Branch(sys+"del_phi", &b_del_phi, sys+"del_phi/F");
    tr->Branch(sys+"del_r", &b_del_r, sys+"del_r/F");
    tr->Branch(sys+"del_r12", &b_del_r12, sys+"del_r12/F");
    tr->Branch(sys+"raw_mass", &b_raw_mass, sys+"raw_mass/F");
    tr->Branch(sys+"jet1_pt", &b_jet1_pt, sys+"jet1_pt/F");
    tr->Branch(sys+"jet1_eta", &b_jet1_eta, sys+"jet1_eta/F");
    tr->Branch(sys+"jet1_phi", &b_jet1_phi, sys+"jet1_phi/F");
    tr->Branch(sys+"jet2_pt", &b_jet2_pt, sys+"jet2_pt/F");
    tr->Branch(sys+"jet2_eta", &b_jet2_eta, sys+"jet2_eta/F");
    tr->Branch(sys+"jet2_phi", &b_jet2_phi, sys+"jet2_phi/F");
    tr->Branch(sys+"jet3_pt", &b_jet3_pt, sys+"jet3_pt/F");
    tr->Branch(sys+"jet3_eta", &b_jet3_eta, sys+"jet3_eta/F");
    tr->Branch(sys+"jet3_phi", &b_jet3_phi, sys+"jet3_phi/F");
    tr->Branch(sys+"met", &b_met, sys+"met/F");
    tr->Branch(sys+"gjet1_pt", &b_gjet1_pt, sys+"gjet1_pt/F");
    tr->Branch(sys+"gjet1_eta", &b_gjet1_eta, sys+"gjet1_eta/F");
    tr->Branch(sys+"gjet1_phi", &b_gjet1_phi, sys+"gjet1_phi/F");
    tr->Branch(sys+"gjet2_pt", &b_gjet2_pt, sys+"gjet2_pt/F");
    tr->Branch(sys+"gjet2_eta", &b_gjet2_eta, sys+"gjet2_eta/F");
    tr->Branch(sys+"gjet2_phi", &b_gjet2_phi, sys+"gjet2_phi/F");
    tr->Branch(sys+"gjet3_pt", &b_gjet3_pt, sys+"gjet3_pt/F");
    tr->Branch(sys+"gjet3_eta", &b_gjet3_eta, sys+"gjet3_eta/F");
    tr->Branch(sys+"gjet3_phi", &b_gjet3_phi, sys+"gjet3_phi/F");
    tr->Branch(sys+"gbeta", &b_gbeta, sys+"gbeta/F");
    tr->Branch(sys+"gdel_eta", &b_gdel_eta, sys+"gdel_eta/F");
    tr->Branch(sys+"gdel_phi", &b_gdel_phi, sys+"gdel_phi/F");
    tr->Branch(sys+"gdel_r", &b_gdel_r, sys+"gdel_r/F");
  }   
  
}
ColorCoherenceAnalyzer::~ColorCoherenceAnalyzer(){}

void ColorCoherenceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  edm::Handle<edm::View<cat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);
  if (!iEvent.getByToken(jetToken_, jets)) return;

  edm::Handle<edm::View<cat::MET> > mets;
  iEvent.getByToken(metToken_, mets);

  edm::Handle<int> vtx;
  iEvent.getByToken(vtxToken_, vtx);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects); 

  vector<cat::Jet> selectedjets = selectJets(jets.product());
  if (selectedjets.size() < 3) return;

  for (unsigned int i = 0; i != sys_e.size(); i++){ 
    resetBr();
    TString sys_name = sys_e[i];

    vector<simpleJet> sortedjets = sortJets(selectedjets, sys_name);
    if (sortedjets.size() < 3 ) return;

    vector<float> beta_results = calBeta(sortedjets); 

    b_nVtx = vtx.product()[0];
    b_nJet = selectedjets.size();
    int hlt_count = 0; 
    if(trigHelper.triggerFired("HLT_Jet40")) {b_hlt_40_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_Jet60")) {b_hlt_60_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_Jet80")) {b_hlt_80_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet140_v")) {b_hlt_140_pass = 1; hlt_count++;}  
    if(trigHelper.triggerFired("HLT_PFJet320_v")) {b_hlt_320_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet400_v")) {b_hlt_400_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet450_v")) {b_hlt_450_pass = 1; hlt_count++;}
    if(trigHelper.triggerFired("HLT_PFJet500_v")) {b_hlt_500_pass = 1; hlt_count++;}
    if (hlt_count < 1) return;
    b_beta = beta_results[0];
    b_del_eta = beta_results[1];
    b_del_phi = beta_results[2];
    b_del_r = beta_results[3];
    b_del_r12 = beta_results[4];
    b_raw_mass = beta_results[5];

    b_jet1_pt = sortedjets[0].pt();
    b_jet1_eta = sortedjets[0].eta();
    b_jet1_phi = sortedjets[0].phi();
    b_jet2_pt = sortedjets[1].pt();
    b_jet2_eta = sortedjets[1].eta();
    b_jet2_phi = sortedjets[1].phi();
    b_jet3_pt = sortedjets[2].pt();
    b_jet3_eta = sortedjets[2].eta();
    b_jet3_phi = sortedjets[2].phi();
    b_met = mets->begin()->et()/mets->begin()->sumEt();
    ttree_[i]->Fill();

  }
}


vector<cat::Jet> ColorCoherenceAnalyzer::selectJets(const edm::View<cat::Jet>* jets)
{
  vector<cat::Jet> seljets;
  for (auto jet : *jets) {
    if (!jet.LooseId()) continue;
    if (jet.pt() <= 20.) continue;
    if (jet.pileupJetId() <0.9) continue;
    seljets.push_back(jet);
  }
  return seljets;
}

simpleJet ColorCoherenceAnalyzer::sysJet(cat::Jet jet, TString sys_name)
{
  simpleJet sJet(jet.pt(), jet.eta(), jet.phi(), jet.p4().M());
  if (sys_name == "eup") sJet.SetPt(jet.pt()*jet.shiftedEnUp());
  if (sys_name == "edown") sJet.SetPt(jet.pt()*jet.shiftedEnDown());
  if (sys_name == "esup") sJet.SetPt(jet.pt()*jet.smearedResUp());
  if (sys_name == "esdown") sJet.SetPt(jet.pt()*jet.smearedResDown());
  if (sys_name == "jar") sJet = jarJet(jet);

  return sJet;
}

simpleJet ColorCoherenceAnalyzer::jarJet(cat::Jet jet)
{
  simpleJet sJet(jet.pt(), jet.eta(), jet.phi(), jet.p4().M());
  
  return sJet;
}

vector<simpleJet> ColorCoherenceAnalyzer::sortJets(vector<cat::Jet> seljets, TString sys_name)
{
  vector<simpleJet> sJet;
  for (auto jet : seljets){
    simpleJet tmp = sysJet(jet, sys_name);
    if (tmp.pt() > 30.0) sJet.push_back(tmp);
  }
  sort(sJet.begin(), sJet.end(), GtByPt());
  return sJet;
}

vector<float> ColorCoherenceAnalyzer::calBeta(vector<simpleJet> seljets)
{
  vector<float> beta_res;
  float del_eta = copysign(1.0, seljets[1].eta())*(seljets[2].eta() - seljets[1].eta());
  float del_phi = reco::deltaPhi(seljets[2].phi(), seljets[1].phi());
  float del_r = reco::deltaR(seljets[2].eta(), seljets[2].phi(), seljets[1].eta(), seljets[1].phi());
  float beta = atan2(del_phi, del_eta);
  float del_r12 = reco::deltaR(seljets[1].eta(), seljets[1].phi(), seljets[0].eta(), seljets[0].phi());
  float raw_mass = (seljets[0] + seljets[1]).M();
  beta_res.push_back(beta);
  beta_res.push_back(del_eta);
  beta_res.push_back(del_phi);
  beta_res.push_back(del_r);
  beta_res.push_back(del_r12);
  beta_res.push_back(raw_mass);

  return beta_res;
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

void ColorCoherenceAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  //  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ColorCoherenceAnalyzer);
