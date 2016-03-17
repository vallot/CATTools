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

#include "TTree.h"
#include "TLorentzVector.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;
using namespace cat;

class h2muAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit h2muAnalyzer(const edm::ParameterSet&);
  ~h2muAnalyzer() {};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  cat::MuonCollection selectMuons(const cat::MuonCollection& muons ) const;
  cat::ElectronCollection selectElecs(const cat::ElectronCollection& elecs ) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, cat::MuonCollection & recolep) const;
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  int preSelect(const cat::JetCollection& seljets, float MET) const;
  int jetCategory(const cat::JetCollection& seljets, float MET, float ll_pt) const;
  int etaCategory(float lep1_eta, float lep2_eta) const;

  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genweightToken_, puweightToken_;
  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  TTree * ttree_;
  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_tri, b_filtered;
  float b_met, b_weight, b_puweight;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;

  int b_cat_jet;
  int b_cat_eta;
  bool b_isLoose, b_isMedium, b_isTight;

  float b_gen_lep1_pt, b_gen_lep1_eta, b_gen_lep1_phi, b_gen_lep1_ptRes;
  bool b_gen_lep1_isLoose, b_gen_lep1_isMedium, b_gen_lep1_isTight;
  float b_gen_lep2_pt, b_gen_lep2_eta, b_gen_lep2_phi, b_gen_lep2_ptRes;
  bool b_gen_lep2_isLoose, b_gen_lep2_isMedium, b_gen_lep2_isTight;
  float b_gen_ll_pt, b_gen_ll_eta, b_gen_ll_phi, b_gen_ll_m;

  bool runOnMC_;
};

h2muAnalyzer::h2muAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  muonToken_ = consumes<cat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("nom", "nom");
  ttree_->Branch("nvertex", &b_nvertex, "nvertex/I");
  ttree_->Branch("step", &b_step, "step/I");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("nbjet", &b_nbjet, "nbjet/I");
  ttree_->Branch("step1", &b_step1, "step1/O");
  ttree_->Branch("step2", &b_step2, "step2/O");
  ttree_->Branch("step3", &b_step3, "step3/O");
  ttree_->Branch("step4", &b_step4, "step4/O");
  ttree_->Branch("step5", &b_step5, "step5/O");
  ttree_->Branch("tri", &b_tri, "tri/O");
  ttree_->Branch("filtered", &b_filtered, "filtered/O");
  ttree_->Branch("met", &b_met, "met/F");
  ttree_->Branch("weight", &b_weight, "weight/F");
  ttree_->Branch("puweight", &b_puweight, "puweight/F");

  ttree_->Branch("lep1_pt", &b_lep1_pt, "lep1_pt/F");
  ttree_->Branch("lep1_eta", &b_lep1_eta, "lep1_eta/F");
  ttree_->Branch("lep1_phi", &b_lep1_phi, "lep1_phi/F");
  ttree_->Branch("lep2_pt", &b_lep2_pt, "lep2_pt/F");
  ttree_->Branch("lep2_eta", &b_lep2_eta, "lep2_eta/F");
  ttree_->Branch("lep2_phi", &b_lep2_phi, "lep2_phi/F");
  ttree_->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
  ttree_->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  ttree_->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  ttree_->Branch("ll_m", &b_ll_m, "ll_m/F");

  ttree_->Branch("isLoose", &b_isLoose, "isLoose/O");
  ttree_->Branch("isMedium", &b_isMedium, "isMedium/O");
  ttree_->Branch("isTight", &b_isTight, "isTight/O");
  
  //final hierachy
  //(e.g. In case of 0,1jet, Tight and Loose.Otherwise 2jet include VBF Tight, ggF Tight, Loose)
  ttree_->Branch("cat_jet", &b_cat_jet, "cat_jet/I");
  //Geometrical Categorization
  //only included 0jet and 1jet
  ttree_->Branch("cat_eta", &b_cat_eta, "cat_eta/I");

  ttree_->Branch("gen_lep1_pt", &b_gen_lep1_pt, "gen_lep1_pt/F");
  ttree_->Branch("gen_lep1_eta", &b_gen_lep1_eta, "gen_lep1_eta/F");
  ttree_->Branch("gen_lep1_phi", &b_gen_lep1_phi, "gen_lep1_phi/F");
  ttree_->Branch("gen_lep1_ptRes", &b_gen_lep1_ptRes, "gen_lep1_ptRes/F");
  ttree_->Branch("gen_lep1_isLoose", &b_gen_lep1_isLoose, "gen_lep1_isLoose/O");
  ttree_->Branch("gen_lep1_isMedium", &b_gen_lep1_isMedium, "gen_lep1_isMedium/O");
  ttree_->Branch("gen_lep1_isTight", &b_gen_lep1_isTight, "gen_lep1_isTight/O");

  ttree_->Branch("gen_lep2_pt", &b_gen_lep2_pt, "gen_lep2_pt/F");
  ttree_->Branch("gen_lep2_eta", &b_gen_lep2_eta, "gen_lep2_eta/F");
  ttree_->Branch("gen_lep2_phi", &b_gen_lep2_phi, "gen_lep2_phi/F");
  ttree_->Branch("gen_lep2_ptRes", &b_gen_lep2_ptRes, "gen_lep2_ptRes/F");
  ttree_->Branch("gen_lep2_isLoose", &b_gen_lep2_isLoose, "gen_lep2_isLoose/O");
  ttree_->Branch("gen_lep2_isMedium", &b_gen_lep2_isMedium, "gen_lep2_isMedium/O");
  ttree_->Branch("gen_lep2_isTight", &b_gen_lep2_isTight, "gen_lep2_isTight/O");

  ttree_->Branch("gen_ll_pt", &b_gen_ll_pt, "gen_ll_pt/F");
  ttree_->Branch("gen_ll_eta", &b_gen_ll_eta, "gen_ll_eta/F");
  ttree_->Branch("gen_ll_phi", &b_gen_ll_phi, "gen_ll_phi/F");
  ttree_->Branch("gen_ll_m", &b_gen_ll_m, "gen_ll_m/F");

}

void h2muAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  b_nvertex = 0;b_step = -1;b_channel = 0;b_njet = 0;b_nbjet = 0;
  b_step1 = 0;b_step2 = 0;b_step3 = 0;b_step4 = 0;b_step5 = 0;b_tri = 0;b_filtered = 0;
  b_met = -9;
  b_weight = 1; b_puweight = 1;

  b_lep1_pt = -9;b_lep1_eta = -9;b_lep1_phi = -9;
  b_lep2_pt = -9;b_lep2_eta = -9;b_lep2_phi = -9;
  b_ll_pt = -9;b_ll_eta = -9;b_ll_phi = -9;b_ll_m = -9;

  b_isLoose = 0; b_isMedium = 0; b_isTight = 0;

  b_cat_jet = -1;
  b_cat_eta = -1;

  b_gen_lep1_pt = 0;b_gen_lep1_eta = 0;b_gen_lep1_phi = 0;b_gen_lep1_ptRes = 0;
  b_gen_lep1_isLoose = 0;b_gen_lep1_isMedium = 0;b_gen_lep1_isTight = 0;
  b_gen_lep2_pt = 0;b_gen_lep2_eta = 0;b_gen_lep2_phi = 0;b_gen_lep2_ptRes = 0;
  b_gen_lep2_isLoose = 0;b_gen_lep2_isMedium = 0;b_gen_lep2_isTight = 0;
  b_gen_ll_pt = 0;b_gen_ll_eta = 0;b_gen_ll_phi = 0;b_gen_ll_m = 0;

  if (runOnMC_){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    b_puweight = *puweightHandle;
    edm::Handle<float> genweightHandle;
    iEvent.getByToken(genweightToken_, genweightHandle);
    b_weight = (*genweightHandle)*b_puweight;
  }
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    ttree_->Fill();
    return;
  }
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  b_nvertex = *nGoodVertexHandle;

  edm::Handle<int> lumiSelectionHandle;
  iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
  if (!runOnMC_){
    if (*lumiSelectionHandle == 0) return;
  }
 
  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;
  // if (!b_filtered){
  //   ttree_->Fill();
  //   return;
  // }
  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);
  
  edm::Handle<reco::GenParticleCollection> genParticles;

  cat::MuonCollection selectedMuons = selectMuons( *muons );
  sort(selectedMuons.begin(), selectedMuons.end(), GtByCandPt());
  
  if (runOnMC_){
    iEvent.getByToken(mcLabel_,genParticles);
    bool bosonSample = false;
    TLorentzVector genLep1;
    TLorentzVector genLep2;
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
        if (g.charge() > 0) genLep1.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), g.mass());
        else genLep2.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), g.mass());
      }
    }
    if (bosonSample){
      b_gen_lep1_pt = genLep1.Pt();b_gen_lep1_eta = genLep1.Eta();b_gen_lep1_phi = genLep1.Phi();
      b_gen_lep2_pt = genLep2.Pt();b_gen_lep2_eta = genLep2.Eta();b_gen_lep2_phi = genLep2.Phi();
      TLorentzVector gen_ll = genLep1 + genLep2;
      b_gen_ll_pt = gen_ll.Pt(); b_gen_ll_eta = gen_ll.Eta(); b_gen_ll_phi = gen_ll.Phi(); b_gen_ll_m = gen_ll.M();

      for (auto m : selectedMuons){
        if (genLep1.DeltaR(m.tlv()) < 0.1){
          b_gen_lep1_ptRes = (m.pt()-genLep1.Pt())/genLep1.Pt();
          b_gen_lep1_isLoose = m.isLooseMuon(); b_gen_lep1_isMedium = m.isMediumMuon(); b_gen_lep1_isTight = m.isTightMuon();
        }
        if (genLep2.DeltaR(m.tlv()) < 0.1){
          b_gen_lep2_ptRes = (m.pt()-genLep2.Pt())/genLep2.Pt();
          b_gen_lep2_isLoose = m.isLooseMuon(); b_gen_lep2_isMedium = m.isMediumMuon(); b_gen_lep2_isTight = m.isTightMuon();
        }
      }
    }
  }
  
  if (selectedMuons.size() < 2){
    ttree_->Fill();
    return;
  }
  //b_muEtaBin = muEtaBin(selectedMuons);
			
  b_step = 1;
  b_step1 = true;
  
  b_lep1_pt = selectedMuons[0].pt(); b_lep1_eta = selectedMuons[0].eta(); b_lep1_phi = selectedMuons[0].phi();
  b_lep2_pt = selectedMuons[1].pt(); b_lep2_eta = selectedMuons[1].eta(); b_lep2_phi = selectedMuons[1].phi();

  b_isLoose = (selectedMuons[0].isLooseMuon() && selectedMuons[1].isLooseMuon());
  b_isMedium = (selectedMuons[0].isMediumMuon() && selectedMuons[1].isMediumMuon());
  b_isTight = (selectedMuons[0].isTightMuon() && selectedMuons[1].isTightMuon());

  TLorentzVector tlv_ll = selectedMuons[0].tlv() + selectedMuons[1].tlv();
  b_ll_pt = tlv_ll.Pt(); b_ll_eta = tlv_ll.Eta(); b_ll_phi = tlv_ll.Phi(); b_ll_m = tlv_ll.M();

  TLorentzVector met = mets->front().tlv();
  b_met = met.Pt();

  cat::JetCollection&& selectedJets = selectJets( *jets.product(), selectedMuons );  
  b_njet = selectedJets.size();
  b_cat_eta = etaCategory(b_lep1_eta, b_lep2_eta);

  int ll_charge = selectedMuons[0].charge()*selectedMuons[1].charge();

  if (ll_charge > 0){
    ttree_->Fill();
    return;
  }
  b_step = 2;
  b_step2 = true;

  if (b_ll_m < 20){
    ttree_->Fill();
    return;
  }
  b_step = 3;
  b_step3 = true;
  
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  if (trigHelper.triggerFired("HLT_IsoMu20_v") || trigHelper.triggerFired("HLT_IsoTrkMu20_v")){
    b_step = 4;
    b_step4 = true;
  }

  if ( trigHelper.triggerMatched("HLT_IsoMu20_v", selectedMuons[0]) ||
       trigHelper.triggerMatched("HLT_IsoMu20_v", selectedMuons[1]) ||
       trigHelper.triggerMatched("HLT_IsoTrkMu20_v", selectedMuons[0]) ||
       trigHelper.triggerMatched("HLT_IsoTrkMu20_v", selectedMuons[1])){
    b_tri = true;
    b_step = 5;
    b_step5 = true;
  }

  // -----------------------------  Jet Category  -----------------------------------
  b_cat_jet = jetCategory(selectedJets, b_met, b_ll_pt);

  ttree_->Fill();
}

cat::MuonCollection h2muAnalyzer::selectMuons(const cat::MuonCollection& muons ) const
{
  cat::MuonCollection selmuons;
  for (auto mu : muons) {
    if (mu.pt() <= 20.) continue;
    if (std::abs(mu.eta()) >= 2.4) continue;
    if (!mu.isLooseMuon()) continue;
    if (mu.relIso(0.4) >= 0.15) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }

  return selmuons;
}

cat::ElectronCollection h2muAnalyzer::selectElecs(const cat::ElectronCollection& elecs ) const
{
  cat::ElectronCollection selelecs;
  for (auto el : elecs) {
    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) continue;
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return selelecs;
}

cat::JetCollection h2muAnalyzer::selectJets(const cat::JetCollection& jets, cat::MuonCollection& recomu ) const
{
  cat::JetCollection seljets;
  for (auto jet : jets) {    
    if (!jet.LooseId()) continue;
    if (jet.pt() <= 30.) continue;
    if (std::abs(jet.eta()) >= 2.4)  continue;
    //if (jet.tlv().DeltaR(recomu[0]) <= 0.4) continue;
    //if (jet.tlv().DeltaR(recomu[1]) <= 0.4) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}

int h2muAnalyzer::preSelect(const cat::JetCollection& seljets, float met) const
{
  int njet = seljets.size();
  if (njet>1){
    if (seljets[0].pt()>40 && seljets[1].pt()>30 && met<40){
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

int h2muAnalyzer::jetCategory(const cat::JetCollection& seljets, float met, float ll_pt) const
{
  int presel = preSelect(seljets, met);
  if (presel==1){
    if (ll_pt<=10){return 1;}
    else{return 2;}
  }
  if (presel==2){
    if (ll_pt<=10){return 3;}
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

int h2muAnalyzer::etaCategory(float lep1_eta, float lep2_eta) const
{
  const float eta_mu[2] = {std::abs(lep1_eta),std::abs(lep2_eta)};
  int GC=0;
  for(int i=0; i<2; i++){
    if      (eta_mu[i] < 0.8) GC += 1;
    else if (eta_mu[i] < 1.5) GC += 10;
    else if (eta_mu[i] < 2.4) GC += 100;
  }
  if (GC == 2)   return 1; // BB
  if (GC == 11)  return 2; // BO
  if (GC == 101) return 3; // BE
  if (GC == 20)  return 4; // OO
  if (GC == 110) return 5; // OE
  if (GC == 200) return 6; // EE
  
  return 0;
}
//define this as a plug-in
DEFINE_FWK_MODULE(h2muAnalyzer);
