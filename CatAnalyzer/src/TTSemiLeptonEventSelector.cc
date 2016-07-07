#include "CATTools/CatAnalyzer/interface/TTSemiLeptonEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;


TTSemiLeptonEventSelector::TTSemiLeptonEventSelector(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iC ): TTEventSelector(iConfig, iC ) {

  
  vetoMuonPtCut_ =  iConfig.getParameter<double>("vetoMuonPtCut");
  vetoMuonEtaCut_ = iConfig.getParameter<double>("vetoMuonEtatCut");
  vetoMuonIsoCut_ = iConfig.getParameter<double>("vetoMuonIsoCut");


  vetoElectronPtCut_  = iConfig.getParameter<double>("vetoElectronPtCut");
  vetoElectronEtaCut_ = iConfig.getParameter<double>("vetoElectronEtatCut");
  vetoElectronIsoCut_ = iConfig.getParameter<double>("vetoElectronIsoCut");
  vetoElectronIDCut_  = iConfig.getParameter<double>("vetoElectronIDCut");

  trigTokenMUJET_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigTokenMUJET"));
  trigTokenELJET_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigTokenMUJET"));


}

void TTSemiLeptonEventSelector::resetBranch(){}

void TTSemiLeptonEventSelector::setBranch(TTree* tree,  int sys ) { 
  return;
}

bool TTSemiLeptonEventSelector::isVetoMuon(cat::Muon mu) const
{
  if (!mu.isLooseMuon() ) return false;
  if (mu.pt() < vetoMuonPtCut_ ) return false;
  if (std::abs(mu.eta())> vetoMuonEtaCut_ ) return false;
  if (mu.relIso(0.4) > vetoMuonIsoCut_ ) return false;
  return true; 
}
bool TTSemiLeptonEventSelector::isSelectMuon(cat::Muon mu) const
{
  if (!mu.isTightMuon() ) return false;
  if (mu.pt() < muonPtCut_ ) return false;
  if (std::abs(mu.eta())> muonEtaCut_ ) return false;
  if (mu.relIso(0.4) > muonIsoCut_ ) return false;
  return true; 
}
bool TTSemiLeptonEventSelector::isVetoElec(cat::Electron el) const
{
  if (el.pt() < vetoElectronPtCut_ ) return false;
  if (std::abs(el.eta()) > vetoElectronEtaCut_) return false;
  if ( !el.electronID(vetoElectronIDCut_) ) return false;
  //if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) return false;

  /*
  if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) return false;
  if ( std::abs(el.scEta()) <= 1.479 ) { if ( el.relIso(0.3) > 0.126) return false; }
  else                                 { if ( el.relIso(0.3) > 0.144) return false; }
  */

  return true;
}

bool TTSemiLeptonEventSelector::isSelectElec(cat::Electron el) const
{
  if (el.pt() < electronPtCut_ ) return false;
  if (std::abs(el.eta()) > electronEtaCut_ ) return false;
  if ( !el.electronID(electronIDCut_) ) return false;

  /*
  if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) return false;
  if ( std::abs(el.scEta()) <= 1.479 ) { if ( el.relIso(0.3) > 0.0766) return false; }
  else                                 { if ( el.relIso(0.3) > 0.0678) return false; }
  */

  return true;
}

float TTSemiLeptonEventSelector::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, cat::MuonCollection& vetomuons, sys_e sys ) const
{
  float weight = 1.;
  for (auto& m : muons) {
    cat::Muon mu(m);
    if (sys == sys_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == sys_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());
 
    if ( isSelectMuon( mu) ) selmuons.push_back(mu);
    else if( isVetoMuon(mu) ) vetomuons.push_back(mu);
  }
  return weight;
}

float TTSemiLeptonEventSelector::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, cat::ElectronCollection& vetoelecs, sys_e sys) const
{
  float weight = 1.;
  for (auto& e : elecs) {
    cat::Electron el(e);
    if (sys == sys_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == sys_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if ( isSelectElec(el)) selelecs.push_back(el);
    else if ( isVetoElec(el)) vetoelecs.push_back(el);
  }
  return weight;
}
 

int TTSemiLeptonEventSelector::eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys){
  const bool runOnMC = !iEvent.isRealData();
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();

  if (sys == sys_nom) evInfo_.cutflow_[0][evInfo_.channel]++;

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
  cat::MuonCollection selMuons , vetoMuons;
  cat::ElectronCollection selElecs, vetoElecs;
  selectMuons(*muons, selMuons, vetoMuons, (sys_e)sys);
  selectElecs(*electrons, selElecs, vetoElecs, (sys_e)sys);
  if ( selMuons.size()+selElecs.size() ==0 ) {
    if (evInfo_.keepTtbarSignal) tree->Fill();
    return -1;
  }
  if (sys == sys_nom) evInfo_.cutflow_[3][evInfo_.channel]++;

  std::vector<const cat::Lepton*> recolep;
  //auto& recolep = evInfo.recolep_;
  recolep.clear();
  for ( const auto& x : selMuons ) recolep.push_back(&x);
  for ( const auto& x : selElecs ) recolep.push_back(&x);

  sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});

  const cat::Lepton& recolep1 = *recolep[0];


  // Determine channel
  const int pdgIdSum = std::abs(recolep1.pdgId()) ;
  if (pdgIdSum == 13) evInfo_.channel = CH_MUJET; // mu+jet
  if (pdgIdSum == 11) evInfo_.channel = CH_ELJET; // el+jet

  evInfo_.mueffweight    = getMuEffSF(recolep1,  0);
  evInfo_.mueffweight_up = getMuEffSF(recolep1, +1);
  evInfo_.mueffweight_dn = getMuEffSF(recolep1, -1);

  evInfo_.eleffweight    = getElEffSF(recolep1,  0);
  evInfo_.eleffweight_up = getElEffSF(recolep1, +1);
  evInfo_.eleffweight_dn = getElEffSF(recolep1, -1);

  // Trigger results
  // No information for Single lepton trigger
  evInfo_.tri = 0;

  edm::Handle<int> trigHandle;
  if      ( evInfo_.channel == CH_ELJET ) iEvent.getByToken(trigTokenELJET_, trigHandle);
  else if ( evInfo_.channel == CH_MUJET ) iEvent.getByToken(trigTokenMUJET_, trigHandle);
  if ( *trigHandle != 0 ){
    //std::cout<<"Active trigger!"<<std::endl;
    evInfo_.tri = 1.0; // Need to Update
    evInfo_.tri_up = 1.0;
    evInfo_.tri_dn = 1.0;
  }
  //else  std::cout<<"No trigger!"<<std::endl;

  evInfo_.lep1 = recolep1.tlv(); evInfo_.lep1_pid = recolep1.pdgId();

  evInfo_.step1 = true;
  evInfo_.step = 1;
  if (sys == sys_nom) evInfo_.cutflow_[4][evInfo_.channel]++;

  if ( evInfo_.channel == CH_MUJET && vetoMuons.size() ==0 ) {
    if ( evInfo_.step1 ) {
      evInfo_.step2 = true;
      evInfo_.step=2 ;
      if (sys == sys_nom) evInfo_.cutflow_[5][evInfo_.channel]++;
    }
  }
  if ( evInfo_.channel == CH_ELJET && vetoElecs.size() ==0 ) {
    if ( evInfo_.step1 ) {
      evInfo_.step2 = true;
      evInfo_.step=2 ;
      if (sys == sys_nom) evInfo_.cutflow_[5][evInfo_.channel]++;
    }
  }


  auto selectedJets  = selectJets(*jets, recolep, (sys_e)sys);
  auto selectedBJets =  selectBJets(selectedJets);

  evInfo_.njet = selectedJets.size();
  evInfo_.nbjet = selectedBJets.size();

  if ( evInfo_.njet >=1 ){
    if (evInfo_.step2){
      evInfo_.step3 = true;
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[6][evInfo_.channel]++;
    }
  }
  if ( evInfo_.njet >=4 ){
    if (evInfo_.step3){
      evInfo_.step4 = true;
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[7][evInfo_.channel]++;
    }
  }
  if ( evInfo_.nbjet >= 1){
    if (evInfo_.step4){
      evInfo_.step5 = true;
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[8][evInfo_.channel]++;
    }
  }
  if ( evInfo_.nbjet >= 2){
    if (evInfo_.step5){
      evInfo_.step6 = true;
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[9][evInfo_.channel]++;
    }
  }
  if ( evInfo_.nbjet == 2){
    if (evInfo_.step6){
      evInfo_.step7 = true;
      ++evInfo_.step;
      if (sys == sys_nom) evInfo_.cutflow_[10][evInfo_.channel]++;
    }
  }

  return 0;
}






