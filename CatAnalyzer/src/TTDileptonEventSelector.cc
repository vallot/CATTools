#include "CATTools/CatAnalyzer/interface/TTDileptonEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;

class TTEventSelector;
TTDileptonEventSelector::TTDileptonEventSelector(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iC ): TTEventSelector(iConfig, iC ) {
  trigTokenMUEL_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));
  muonIDCut_ = iConfig.getParameter<string>("MuonIDCut");
}

void TTDileptonEventSelector::setBranch(TTree* tree,  int sys ) { 
  return;
}

void TTDileptonEventSelector::resetBranch(){}
 
int TTDileptonEventSelector::eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys){
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

  if (selectedJets.size() >1 ){
    evInfo_.step3 = true;
    if (evInfo_.step == 2){
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[6][evInfo_.channel]++;
    }
  }

  if ((evInfo_.channel == CH_MUEL) || (evInfo_.met > 40.)){
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


