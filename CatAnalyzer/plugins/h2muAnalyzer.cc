#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
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

#include "CATTools/CatAnalyzer/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

class h2muAnalyzer : public edm::EDAnalyzer {
public:
  explicit h2muAnalyzer(const edm::ParameterSet&);
  ~h2muAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
  vector<cat::Muon> selectMuons(const edm::View<cat::Muon>* muons );
  vector<cat::Electron> selectElecs(const edm::View<cat::Electron>* elecs );
  vector<cat::Jet> selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep);
  vector<cat::Jet> selectBJets(vector<cat::Jet> & jets );
  float passingSteps(int channel, float met, float ll_mass, float ll_charge, int selectedJets_size);
  int preSelect(vector<cat::Jet> seljets, float MET);
  int JetCategory(vector<cat::Jet> seljets, float MET, float ll_pt);
  int JetCat_GC(float mu1_eta, float mu2_eta);

  edm::EDGetTokenT<edm::View<cat::Muon> >     muonToken_;
  edm::EDGetTokenT<edm::View<cat::Electron> > elecToken_;
  edm::EDGetTokenT<edm::View<cat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> >      metToken_;
  edm::EDGetTokenT<reco::VertexCollection >   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  TTree * ttree_;
  TTree * t2;

  int b_njet, b_step, b_channel;
  float b_MET;
  float b_lep1_pt, b_lep1_eta, b_lep1_phi;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;
  int b_jetcat_f_hier;  
  int b_jetcat_GC;
  bool b_isMedium, b_isTight;

  float b_gen_lep_pt, b_gen_lep_eta, b_gen_lep_phi;
  float b_reco_lep_pt, b_reco_lep_eta, b_reco_lep_phi;
  float b_resolution;  
  bool b_lep_isLoose, b_lep_isMedium, b_lep_isTight;

  bool runOnMC_;
};
h2muAnalyzer::h2muAnalyzer(const edm::ParameterSet& iConfig)
{
  muonToken_ = consumes<edm::View<cat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<edm::View<cat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<edm::View<cat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<edm::View<cat::MET> >(iConfig.getParameter<edm::InputTag>("mets"));     
  vtxToken_  = consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("tree", "tree");
  t2 = fs->make<TTree>("tree2","tree2");
  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("step", &b_step, "step/I");

  ttree_->Branch("isMedium", &b_isMedium, "isMedium/B");
  ttree_->Branch("isTight", &b_isTight, "isTight/B");

  ttree_->Branch("mu1_pt", &b_lep1_pt, "mu1_pt/F");
  ttree_->Branch("mu1_eta", &b_lep1_eta, "mu1_eta/F");
  ttree_->Branch("mu1_phi", &b_lep1_phi, "mu1_phi/F");

  ttree_->Branch("mu2_pt", &b_lep2_pt, "mu2_pt/F");
  ttree_->Branch("mu2_eta", &b_lep2_eta, "mu2_eta/F");
  ttree_->Branch("mu2_phi", &b_lep2_phi, "mu2_phi/F");

  ttree_->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
  ttree_->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  ttree_->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  ttree_->Branch("ll_m", &b_ll_m, "ll_m/F");

  //final hierachy
  //(e.g. In case of 0,1jet, Tight and Loose.Otherwise 2jet include VBF Tight, ggF Tight, Loose)
  ttree_->Branch("jetcat_f_hier", &b_jetcat_f_hier, "jetcat_f_hier/I");
  
  //Geometrical Categorization
  //only included 0jet and 1jet
  ttree_->Branch("jetcat_GC", &b_jetcat_GC, "jetcat_GC/I");

  //tree2. we should use this for efficiency of recomuon per genmuon.
  t2 ->Branch("gen_lep_pt", &b_gen_lep_pt, "gen_lep_pt/F");
  t2 ->Branch("gen_lep_eta", &b_gen_lep_eta, "gen_lep_eta/F");
  t2 ->Branch("gen_lep_phi", &b_gen_lep_phi, "gen_lep_phi/F");
  t2 ->Branch("reco_lep_pt", &b_reco_lep_pt, "reco_lep_pt/F");
  t2 ->Branch("reco_lep_eta", &b_reco_lep_eta, "reco_lep_eta/F");
  t2 ->Branch("reco_lep_phi", &b_reco_lep_phi, "reco_lep_phi/F");
  t2 ->Branch("resolution", &b_resolution, "resolution/F");
  t2 ->Branch("lep_isLoose", &b_lep_isLoose, "lep_isLoose/B");
  t2 ->Branch("lep_isMedium", &b_lep_isMedium, "lep_isMedium/B");
  t2 ->Branch("lep_isTight", &b_lep_isTight, "lep_isTight/B");
 
}
h2muAnalyzer::~h2muAnalyzer(){}

void h2muAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle<edm::View<cat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<edm::View<cat::Electron> > electrons;
  iEvent.getByToken(elecToken_, electrons);

  edm::Handle<edm::View<cat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<edm::View<cat::MET> > mets;
  iEvent.getByToken(metToken_, mets);
 
  edm::Handle<reco::GenParticleCollection> genParticles;

  b_njet = -1; b_step = 0; b_channel = -1;
  b_MET = -1;
  b_lep1_pt = -9; b_lep1_eta = -9; b_lep1_phi = -9;
  b_lep2_pt = -9; b_lep2_eta = -9; b_lep2_phi = -9;
  b_ll_pt = -9; b_ll_eta = -9; b_ll_phi = -9; b_ll_m = -9;
  b_isMedium = 0; b_isTight = 0;

  b_jetcat_f_hier = -9;
  b_jetcat_GC = -9;  

  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );

  if (runOnMC_){
    iEvent.getByToken(mcLabel_,genParticles);
    for (const reco::GenParticle & g : *genParticles){
      bool isfromZboson = false;
      if (abs(g.pdgId())!=13 || g.pt()<=20.){
        continue;
      }
      for (unsigned int i = 0; i < g.numberOfMothers(); ++i){
	if (g.mother(i)->pdgId()  == 23){ //In case of pdgId() = 23, indicate Z-boson. if it's 25, that becomes higgs.
	  isfromZboson = true;
	}
      }
      
      b_gen_lep_pt = -9; b_gen_lep_eta = -9; b_gen_lep_phi = -9;
      b_reco_lep_pt = -9; b_reco_lep_eta = -9; b_reco_lep_phi = -9;
      b_resolution = -9;  
      b_lep_isLoose = 0; b_lep_isMedium = 0; b_lep_isTight = 0;  
 
      if (!isfromZboson) {
        t2->Fill();
        continue;
      }
      b_gen_lep_pt = g.pt();
      b_gen_lep_eta = g.eta();
      b_gen_lep_phi = g.phi();
      
      for (auto m : selectedMuons){
        b_lep_isLoose = m.isLooseMuon();
        b_lep_isMedium = m.isMediumMuon();
        b_lep_isTight = m.isTightMuon();
        float dr = reco::deltaR(g.eta(), g.phi(), m.eta(), m.phi());
        if (dr < 0.1){
          b_reco_lep_pt = g.pt();
          b_reco_lep_eta = g.eta();
          b_reco_lep_phi = g.phi();
          b_resolution = (m.pt()-g.pt())/g.pt();
          break;
        }
      }
      t2->Fill();
    }    
  }

  if (selectedMuons.size() < 2){
    ttree_->Fill();
    return;
  }

  b_lep1_pt = selectedMuons[0].pt();
  b_lep1_eta = selectedMuons[0].eta();
  b_lep1_phi = selectedMuons[0].phi();

  b_lep2_pt = selectedMuons[1].pt();
  b_lep2_eta = selectedMuons[1].eta();
  b_lep2_phi = selectedMuons[1].phi();

  b_isMedium = (selectedMuons[0].isMediumMuon() && selectedMuons[1].isMediumMuon());
  b_isTight = (selectedMuons[0].isTightMuon() && selectedMuons[1].isTightMuon());
  

  TLorentzVector tlv_ll = selectedMuons[0].tlv() + selectedMuons[1].tlv();
  
  b_ll_pt = tlv_ll.Pt();
  b_ll_eta = tlv_ll.Eta();
  b_ll_phi = tlv_ll.Phi();
  b_ll_m = tlv_ll.M();

  b_step = 1;
  int ll_charge = selectedMuons[0].charge()*selectedMuons[1].charge();

  if (ll_charge < 0)
    b_step = 2;

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);  
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  //  if (triggerFired(triggerNames, triggerBits, "HLT_IsoMu24_eta2p1_v") ){
  if (trigHelper.triggerFired("HLT_IsoMu24_eta2p1_v") ){
    b_step = 3;
    cout << "trigger fired"<<endl;
  }
  // if ( triggerMatched(triggerNames, triggerObjects, "HLT_IsoMu24_eta2p1_v", selectedMuons[0] )
  //      || triggerMatched(triggerNames, triggerObjects, "HLT_IsoMu24_eta2p1_v", selectedMuons[1] ))
  if ( trigHelper.triggerMatched("HLT_IsoMu24_eta2p1_v", selectedMuons[0] )
       || trigHelper.triggerMatched("HLT_IsoMu24_eta2p1_v", selectedMuons[1] ))
    b_step = 4;

  TLorentzVector met = mets->front().tlv();
  b_MET = met.Pt();

  vector<TLorentzVector> recolep; 
  //for (auto lep : selectedMuons){ recolep.push_back(lep.tlv()); }
  //vector<cat::Electron> selectedElectrons = selectElecs( electrons.product() );
  vector<cat::Jet> selectedJets = selectJets( jets.product(), recolep );

  b_njet = selectedJets.size();

  //  float step = passingSteps( channel, met.Pt(), (recolep[0]+recolep[1]).M(), ll_charge, selectedJets.size() );
  
  // -----------------------------  Jet Category  -----------------------------------
  b_jetcat_f_hier = JetCategory(selectedJets, b_MET, b_ll_pt);
  b_jetcat_GC = JetCat_GC(b_lep1_eta, b_lep2_eta);
   
  ttree_->Fill();
}

vector<cat::Muon> h2muAnalyzer::selectMuons(const edm::View<cat::Muon>* muons )
{
  vector<cat::Muon> selmuons;
  for (auto mu : *muons) {
    if (!mu.isLooseMuon()) continue;
    if (mu.pt() <= 20.) continue;
    if (fabs(mu.eta()) >= 2.4) continue;
    if (mu.relIso(0.4) >= 0.12) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }

  return selmuons;
}

vector<cat::Electron> h2muAnalyzer::selectElecs(const edm::View<cat::Electron>* elecs )
{
  vector<cat::Electron> selelecs;
  for (auto el : *elecs) {
    if (!el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium")) continue;
    if (!el.passConversionVeto()) continue;
    if (!el.isPF()) continue;
    if (el.pt() <= 20.) continue;
    if ((fabs(el.scEta()) <= 1.4442) && (el.relIso(0.3) >= 0.1649)) continue;
    if ((fabs(el.scEta()) >= 1.566) && (el.relIso(0.3) >= 0.2075)) continue;
    if ((fabs(el.scEta()) > 1.4442) && (fabs(el.scEta()) < 1.566)) continue;
    if (fabs(el.eta()) >= 2.5) continue;
    if (el.pt() < 5) continue;
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return selelecs;
}

vector<cat::Jet> h2muAnalyzer::selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep )
{
  vector<cat::Jet> seljets;
  for (auto jet : *jets) {
    if (!jet.LooseId()) continue;
    if (jet.pt() <= 30.) continue;
    if (fabs(jet.eta()) >= 2.4)	continue;
    //if (jet.tlv().DeltaR(recolep[0]) <= 0.4) continue;
    //if (jet.tlv().DeltaR(recolep[1]) <= 0.4) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}


float h2muAnalyzer::passingSteps(int channel, float met, float ll_mass, float ll_charge, int selectedJets_size)
{
  int step = 0;
  if (ll_mass <= 20.) return step;
  if (ll_charge > 0.) return step;
  step = 1;
  if (channel != 1){
    if ((ll_mass > 76) and (ll_mass < 106)) return step;
  }
  step = 2;
  if (selectedJets_size < 2) return step;
  step = 3;
  if (channel == 1){
    step = 4;
  }
  else{
    if (met <= 40.) return step;
  }
  step = 4;

  return step;
}

int h2muAnalyzer::preSelect(vector<cat::Jet> seljets, float MET)
{
  int njet = seljets.size();
  if (njet>1){
    if (seljets[0].pt()>40 && seljets[1].pt()>30 && MET<40){
      return 3;
    }    
  }
  if (njet==1){
    return 2;
  }
  if (njet==0){
    return 1;
  }
  return 0;
}

int h2muAnalyzer::JetCategory(vector<cat::Jet> seljets, float MET, float ll_pt)
{
  int presel = preSelect(seljets, MET);
  if (presel==1){
    if (b_ll_pt<=10){return 1;}
    else{return 2;}
  }
  if (presel==2){
    if (b_ll_pt<=10){return 3;}
    else{return 4;}
  }
  if (presel==3){
    TLorentzVector M_jets = seljets[0].tlv() + seljets[1].tlv();
    auto delta_eta = seljets[0].eta()-seljets[1].eta();
    bool VBF_Tight = (M_jets.M() > 650 && abs(delta_eta) > 3.5);
    bool ggF_Tight = (M_jets.M() > 250 && ll_pt > 50);
    if (VBF_Tight || ggF_Tight){
      if (!ggF_Tight){return 5;} //not ggF_Tight but only VBF_Tight
      if (!VBF_Tight){return 6;}//also contrast of above
      if (VBF_Tight && ggF_Tight){return 7;}
    }
    else {return 8;}
  }
  return 0;
}

int h2muAnalyzer::JetCat_GC(float mu1_eta, float mu2_eta)
{
  float eta_mu[2] = {abs(mu1_eta),abs(mu2_eta)};
  float GC=0;
  for(int i=0; i<2; i++){
    if (eta_mu[i] < 0.8) GC += 1;
    if (eta_mu[i] > 0.8 && eta_mu[i] < 1.5) GC += 10;
    if (eta_mu[i] > 1.5 && eta_mu[i] < 2.4) GC += 100;
  }  
  return GC;
}

void h2muAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(h2muAnalyzer);
