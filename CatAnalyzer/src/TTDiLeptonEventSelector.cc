#include "CATTools/CatAnalyzer/interface/TTDiLeptonEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;

class TTEventSelector;
TTDiLeptonEventSelector::TTDiLeptonEventSelector(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iC ): TTEventSelector(iConfig, iC ) {
  std::cout<<"Running for TTbar dilepton decay channel Reference Event Selection!"<<std::endl;
  trigTokenMUEL_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));
  muonIDCut_ = iConfig.getParameter<string>("MuonIDCut");
  
}

void TTDiLeptonEventSelector::setBranch(TTree* tree,  int sys ) { 
  return;
}

void TTDiLeptonEventSelector::resetBranch(){}
 
int TTDiLeptonEventSelector::eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys){
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  const bool runOnMC = !iEvent.isRealData();
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    if (evInfo_.keepTtbarSignal) tree->Fill();
    std::cout<<"No PV"<<std::endl;
    return -1;
  }
  if (sys == sys_nom) evInfo_.cutflow_[1][evInfo_.channel]++;

  // const reco::Vertex &PV = vertices->front();
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  evInfo_.nvertex = *nGoodVertexHandle;

  edm::Handle<int> lumiSelectionHandle;
  iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
  if (!runOnMC){
    if (*lumiSelectionHandle == 0) return -2;  // Critical problem! Terminate the job for this sample.
  }

  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  evInfo_.filtered = *recoFiltersHandle == 0 ? false : true;
  if (sys == sys_nom) evInfo_.cutflow_[2][evInfo_.channel]++;

  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);

  // Find leptons and sort by pT
  cat::MuonCollection selMuons;
  cat::ElectronCollection selElecs;
  selectMuons(*muons, selMuons, (sys_e)sys);
  selectElecs(*electrons, selElecs, (sys_e)sys);
  if ( selMuons.size()+selElecs.size() < 2 ) {
    if (evInfo_.keepTtbarSignal) tree->Fill();
    return -1;
  }
  if (sys == sys_nom) evInfo_.cutflow_[3][evInfo_.channel]++;

  //std::vector<const cat::Lepton*> recolep;
  auto& recolep = evInfo_.recolep_;
  recolep.clear();
  for ( const auto& x : selMuons ) recolep.push_back(&x);
  for ( const auto& x : selElecs ) recolep.push_back(&x);

  sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
  //evInfo_.is3lep = recolep.size();
  recolep.erase(recolep.begin()+2,recolep.end());
  if ( recolep[0]->charge() < 0 ){
    swap(recolep[0], recolep[1]);
  }
  const cat::Lepton& recolep1 = *recolep[0];
  const cat::Lepton& recolep2 = *recolep[1];

  // Determine channel
  const int pdgIdSum = std::abs(recolep1.pdgId()) + std::abs(recolep2.pdgId());
  if (pdgIdSum == 24) evInfo_.channel = CH_MUEL; // emu
  if (pdgIdSum == 22) evInfo_.channel = CH_ELEL; // ee
  if (pdgIdSum == 26) evInfo_.channel = CH_MUMU; // mumu

  evInfo_.mueffweight    = getMuEffSF(recolep1,  0)*getMuEffSF(recolep2,  0);
  evInfo_.mueffweight_up = getMuEffSF(recolep1, +1)*getMuEffSF(recolep2, +1);
  evInfo_.mueffweight_dn = getMuEffSF(recolep1, -1)*getMuEffSF(recolep2, -1);

  evInfo_.eleffweight    = getElEffSF(recolep1,  0)*getElEffSF(recolep2,  0);
  evInfo_.eleffweight_up = getElEffSF(recolep1, +1)*getElEffSF(recolep2, +1);
  evInfo_.eleffweight_dn = getElEffSF(recolep1, -1)*getElEffSF(recolep2, -1);

  // Trigger results
  edm::Handle<int> trigHandle;
  if      ( evInfo_.channel == CH_ELEL ) iEvent.getByToken(trigTokenELEL_, trigHandle);
  else if ( evInfo_.channel == CH_MUMU ) iEvent.getByToken(trigTokenMUMU_, trigHandle);
  else if ( evInfo_.channel == CH_MUEL ) iEvent.getByToken(trigTokenMUEL_, trigHandle);
  if ( *trigHandle != 0 ) {
    evInfo_.tri = computeTrigSF(recolep1, recolep2);
    evInfo_.tri_up = computeTrigSF(recolep1, recolep2,  1);
    evInfo_.tri_dn = computeTrigSF(recolep1, recolep2, -1);
  }

  evInfo_.lep1 = recolep1.tlv(); evInfo_.lep1_pid = recolep1.pdgId();
  evInfo_.lep2 = recolep2.tlv(); evInfo_.lep2_pid = recolep2.pdgId();
  evInfo_.dilep = evInfo_.lep1+evInfo_.lep2;
  const auto tlv_ll = recolep1.p4()+recolep2.p4();

  if (tlv_ll.M() < 20. || recolep1.charge() * recolep2.charge() > 0){
    if (evInfo_.keepTtbarSignal) tree->Fill();
    return -1;
  }
  evInfo_.step1 = true;
  evInfo_.step = 1;
  if (sys == sys_nom) evInfo_.cutflow_[4][evInfo_.channel]++;

  if ( (evInfo_.channel == CH_MUEL) || ((tlv_ll.M() < 76) || (tlv_ll.M() > 106)) ){
    evInfo_.step2 = true;
    evInfo_.step = 2;
    if (sys == sys_nom) evInfo_.cutflow_[5][evInfo_.channel]++;
  }

  //selectedJets.clear();
  //selectedBJets.clear();

  auto selectedJets  = selectJets(*jets, recolep, (sys_e)sys);
  auto selectedBJets =  selectBJets(selectedJets);

  //met.clear();
  evInfo_.metlv  = mets->front().p4();
  evInfo_.met = evInfo_.metlv.pt();
  evInfo_.njet = selectedJets.size();
  evInfo_.nbjet = evInfo_.selectedBJets.size();


  if ((evInfo_.channel == CH_MUEL) || (evInfo_.met > 40.)){
    evInfo_.step3 = true;
    if (evInfo_.step == 2){
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[6][evInfo_.channel]++;
    }
  }
  if (selectedJets.size() >1 ){
    evInfo_.step4 = true;
    if (evInfo_.step == 3){
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[7][evInfo_.channel]++;
    }
  }

  if (selectedBJets.size() > 0){
    evInfo_.step5 = true;
    if (evInfo_.step == 4){
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[8][evInfo_.channel]++;
    }
  }
  return 0;
}
float TTDiLeptonEventSelector::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, sys_e sys) const
{
  float weight = 1.;
  for (auto& m : muons) {
    cat::Muon mu(m);
    if (sys == sys_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == sys_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());

    if (mu.pt() < muonPtCut_) continue;
    if (std::abs(mu.eta()) > std::abs(muonEtaCut_) ) continue;
    // ID Cut for muon. muIDCut_ mus be lower characters.
    if ( muonIDCut_.compare(string("tight"))==0 && !mu.isTightMuon()   ) continue;  
    else if ( muonIDCut_.compare(string("loose"))==0 && !mu.isLooseMuon()   ) continue;  

    if (mu.relIso(0.4) > muonIsoCut_ ) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    //weight *= mu.scaleFactor("NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1");
    //weight *= mu.scaleFactor("NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1");
    selmuons.push_back(mu);
  }
  return weight;
}

float TTDiLeptonEventSelector::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
{
  float weight = 1.;
  for (auto& e : elecs) {
    cat::Electron el(e);
    if (sys == sys_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == sys_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if (el.pt() < electronPtCut_) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > electronEtaCut_ ) continue;
    if ( !el.electronID(electronIDCut_) ) continue;
    //if ( !el.isTrigMVAValid() or !el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") ) continue;
    if (el.relIso(0.3) > electronIsoCut_ ) continue;

    //weight *= el.scaleFactor("mvaEleID-Spring15-25ns-Trig-V1-wp90");
    //weight *= el.scaleFactor("cutBasedElectronID-Spring15-25ns-V1-standalone-medium");
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return weight;
}

