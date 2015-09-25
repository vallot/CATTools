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

  int b_njet, b_step;
  float b_MET;
  float b_mu1_pt, b_mu1_eta, b_mu1_phi;
  float b_mu2_pt, b_mu2_eta, b_mu2_phi;
  float b_diMu_pt, b_diMu_eta, b_diMu_phi, b_diMu_m;
  int b_jetcat_f_hier;  
  int b_jetcat_GC;
  bool b_isLoose, b_isMedium, b_isTight;

  float b_gen_mu1_pt, b_gen_mu1_eta, b_gen_mu1_phi, b_gen_mu1_ptRes;
  bool b_gen_mu1_isLoose, b_gen_mu1_isMedium, b_gen_mu1_isTight;
  float b_gen_mu2_pt, b_gen_mu2_eta, b_gen_mu2_phi, b_gen_mu2_ptRes;
  bool b_gen_mu2_isLoose, b_gen_mu2_isMedium, b_gen_mu2_isTight;
  float b_gen_diMu_pt, b_gen_diMu_eta, b_gen_diMu_phi, b_gen_diMu_m;

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
  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("step", &b_step, "step/I");

  ttree_->Branch("isLoose", &b_isLoose, "isLoose/B");
  ttree_->Branch("isMedium", &b_isMedium, "isMedium/B");
  ttree_->Branch("isTight", &b_isTight, "isTight/B");

  ttree_->Branch("mu1_pt", &b_mu1_pt, "mu1_pt/F");
  ttree_->Branch("mu1_eta", &b_mu1_eta, "mu1_eta/F");
  ttree_->Branch("mu1_phi", &b_mu1_phi, "mu1_phi/F");

  ttree_->Branch("mu2_pt", &b_mu2_pt, "mu2_pt/F");
  ttree_->Branch("mu2_eta", &b_mu2_eta, "mu2_eta/F");
  ttree_->Branch("mu2_phi", &b_mu2_phi, "mu2_phi/F");

  ttree_->Branch("diMu_pt", &b_diMu_pt, "diMu_pt/F");
  ttree_->Branch("diMu_eta", &b_diMu_eta, "diMu_eta/F");
  ttree_->Branch("diMu_phi", &b_diMu_phi, "diMu_phi/F");
  ttree_->Branch("diMu_m", &b_diMu_m, "diMu_m/F");

  //final hierachy
  //(e.g. In case of 0,1jet, Tight and Loose.Otherwise 2jet include VBF Tight, ggF Tight, Loose)
  ttree_->Branch("jetcat_f_hier", &b_jetcat_f_hier, "jetcat_f_hier/I");
  //Geometrical Categorization
  //only included 0jet and 1jet
  ttree_->Branch("jetcat_GC", &b_jetcat_GC, "jetcat_GC/I");

  ttree_->Branch("gen_mu1_pt", &b_gen_mu1_pt, "gen_mu1_pt/F");
  ttree_->Branch("gen_mu1_eta", &b_gen_mu1_eta, "gen_mu1_eta/F");
  ttree_->Branch("gen_mu1_phi", &b_gen_mu1_phi, "gen_mu1_phi/F");
  ttree_->Branch("gen_mu1_ptRes", &b_gen_mu1_ptRes, "gen_mu1_ptRes/F");
  ttree_->Branch("gen_mu1_isLoose", &b_gen_mu1_isLoose, "gen_mu1_isLoose/F");
  ttree_->Branch("gen_mu1_isMedium", &b_gen_mu1_isMedium, "gen_mu1_isMedium/F");
  ttree_->Branch("gen_mu1_isTight", &b_gen_mu1_isTight, "gen_mu1_isTight/F");

  ttree_->Branch("gen_mu2_pt", &b_gen_mu2_pt, "gen_mu2_pt/F");
  ttree_->Branch("gen_mu2_eta", &b_gen_mu2_eta, "gen_mu2_eta/F");
  ttree_->Branch("gen_mu2_phi", &b_gen_mu2_phi, "gen_mu2_phi/F");
  ttree_->Branch("gen_mu2_ptRes", &b_gen_mu2_ptRes, "gen_mu2_ptRes/F");
  ttree_->Branch("gen_mu2_isLoose", &b_gen_mu2_isLoose, "gen_mu2_isLoose/F");
  ttree_->Branch("gen_mu2_isMedium", &b_gen_mu2_isMedium, "gen_mu2_isMedium/F");
  ttree_->Branch("gen_mu2_isTight", &b_gen_mu2_isTight, "gen_mu2_isTight/F");

  ttree_->Branch("gen_diMu_pt", &b_gen_diMu_pt, "gen_diMu_pt/F");
  ttree_->Branch("gen_diMu_eta", &b_gen_diMu_eta, "gen_diMu_eta/F");
  ttree_->Branch("gen_diMu_phi", &b_gen_diMu_phi, "gen_diMu_phi/F");
  ttree_->Branch("gen_diMu_m", &b_gen_diMu_m, "gen_diMu_m/F");
  
}
h2muAnalyzer::~h2muAnalyzer(){}

void h2muAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  b_njet = 0; b_step = 0; b_MET = -9;
  b_mu1_pt = -9; b_mu1_eta = -9; b_mu1_phi = -9;
  b_mu2_pt = -9; b_mu2_eta = -9; b_mu2_phi = -9;
  b_diMu_pt = -9; b_diMu_eta = -9; b_diMu_phi = -9; b_diMu_m = -9;
  b_isLoose = 0; b_isMedium = 0; b_isTight = 0;

  b_jetcat_f_hier = 0;
  b_jetcat_GC = 0;

  b_gen_mu1_pt = 0;b_gen_mu1_eta = 0;b_gen_mu1_phi = 0;b_gen_mu1_ptRes = 0;
  b_gen_mu1_isLoose = 0;b_gen_mu1_isMedium = 0;b_gen_mu1_isTight = 0;
  b_gen_mu2_pt = 0;b_gen_mu2_eta = 0;b_gen_mu2_phi = 0;b_gen_mu2_ptRes = 0;
  b_gen_mu2_isLoose = 0;b_gen_mu2_isMedium = 0;b_gen_mu2_isTight = 0;
  b_gen_diMu_pt = 0;b_gen_diMu_eta = 0;b_gen_diMu_phi = 0;b_gen_diMu_m = 0;
  
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
  
  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );

  if (runOnMC_){
    iEvent.getByToken(mcLabel_,genParticles);
    bool bosonSample = false;
    TLorentzVector genMu1;
    TLorentzVector genMu2;
    for (const reco::GenParticle & g : *genParticles){
      if (abs(g.pdgId())!=13){
        continue;
      }
      bool isfromBoson = false;
      for (unsigned int i = 0; i < g.numberOfMothers(); ++i){
	//In case of pdgId() = 23, indicate Z-boson. if it's 25, that becomes higgs.
	if (g.mother(i)->pdgId() == 23 || g.mother(i)->pdgId() == 25){
	  isfromBoson = true;
	  bosonSample = true;
	}
      }
      if (isfromBoson){
	if (g.charge() > 0) genMu1.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), g.mass());
	else genMu2.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), g.mass());
      }
    }
    if (bosonSample){
      b_gen_mu1_pt = genMu1.Pt();b_gen_mu1_eta = genMu1.Eta();b_gen_mu1_phi = genMu1.Phi();
      b_gen_mu2_pt = genMu2.Pt();b_gen_mu2_eta = genMu2.Eta();b_gen_mu2_phi = genMu2.Phi();
      TLorentzVector gen_diMu = genMu1 + genMu2;
      b_gen_diMu_pt = gen_diMu.Pt(); b_gen_diMu_eta = gen_diMu.Eta(); b_gen_diMu_phi = gen_diMu.Phi(); b_gen_diMu_m = gen_diMu.M();

      for (auto m : selectedMuons){
	if (genMu1.DeltaR(m.tlv()) < 0.1){
	  b_gen_mu1_ptRes = (m.pt()-genMu1.Pt())/genMu1.Pt();      
	  b_gen_mu1_isLoose = m.isLooseMuon(); b_gen_mu1_isMedium = m.isMediumMuon(); b_gen_mu1_isTight = m.isTightMuon();
	}
	if (genMu2.DeltaR(m.tlv()) < 0.1){
	  b_gen_mu2_ptRes = (m.pt()-genMu2.Pt())/genMu2.Pt();      
	  b_gen_mu2_isLoose = m.isLooseMuon(); b_gen_mu2_isMedium = m.isMediumMuon(); b_gen_mu2_isTight = m.isTightMuon();
	}
      }
    }
  }

  if (selectedMuons.size() < 2){
    ttree_->Fill();
    return;
  }
  b_step = 1;

  b_mu1_pt = selectedMuons[0].pt(); b_mu1_eta = selectedMuons[0].eta(); b_mu1_phi = selectedMuons[0].phi();
  b_mu2_pt = selectedMuons[1].pt(); b_mu2_eta = selectedMuons[1].eta(); b_mu2_phi = selectedMuons[1].phi();

  b_isLoose = (selectedMuons[0].isLooseMuon() && selectedMuons[1].isLooseMuon());
  b_isMedium = (selectedMuons[0].isMediumMuon() && selectedMuons[1].isMediumMuon());
  b_isTight = (selectedMuons[0].isTightMuon() && selectedMuons[1].isTightMuon());
  
  TLorentzVector tlv_ll = selectedMuons[0].tlv() + selectedMuons[1].tlv();
  b_diMu_pt = tlv_ll.Pt(); b_diMu_eta = tlv_ll.Eta(); b_diMu_phi = tlv_ll.Phi(); b_diMu_m = tlv_ll.M();

  TLorentzVector met = mets->front().tlv();
  b_MET = met.Pt();

  vector<TLorentzVector> recomu; 
  vector<cat::Jet> selectedJets = selectJets( jets.product(), recomu );

  b_njet = selectedJets.size();
  
  int diMu_charge = selectedMuons[0].charge()*selectedMuons[1].charge();

  if (diMu_charge > 0){
    ttree_->Fill();
    return;
  }
  b_step = 2;
    
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);  
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  if (!trigHelper.triggerFired("HLT_IsoMu24_eta2p1_v")){
  }
  b_step = 3;
  
  if ( !trigHelper.triggerMatched("HLT_IsoMu24_eta2p1_v", selectedMuons[0] )
       && !trigHelper.triggerMatched("HLT_IsoMu24_eta2p1_v", selectedMuons[1] )){
    ttree_->Fill();
    return;
  }
  b_step = 4;
  
  // -----------------------------  Jet Category  -----------------------------------
  b_jetcat_f_hier = JetCategory(selectedJets, b_MET, b_diMu_pt);
  b_jetcat_GC = JetCat_GC(b_mu1_eta, b_mu2_eta);
   
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

vector<cat::Jet> h2muAnalyzer::selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recomu )
{
  vector<cat::Jet> seljets;
  for (auto jet : *jets) {
    if (!jet.LooseId()) continue;
    if (jet.pt() <= 30.) continue;
    if (fabs(jet.eta()) >= 2.4)	continue;
    //if (jet.tlv().DeltaR(recomu[0]) <= 0.4) continue;
    //if (jet.tlv().DeltaR(recomu[1]) <= 0.4) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
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

int h2muAnalyzer::JetCategory(vector<cat::Jet> seljets, float MET, float diMu_pt)
{
  int presel = preSelect(seljets, MET);
  if (presel==1){
    if (b_diMu_pt<=10){return 1;}
    else{return 2;}
  }
  if (presel==2){
    if (b_diMu_pt<=10){return 3;}
    else{return 4;}
  }
  if (presel==3){
    TLorentzVector M_jets = seljets[0].tlv() + seljets[1].tlv();
    auto delta_eta = seljets[0].eta()-seljets[1].eta();
    bool VBF_Tight = (M_jets.M() > 650 && abs(delta_eta) > 3.5);
    bool ggF_Tight = (M_jets.M() > 250 && diMu_pt > 50);
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
