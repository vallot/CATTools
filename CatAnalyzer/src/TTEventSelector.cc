#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;


TTEventSelector::TTEventSelector(const edm::ParameterSet& iConfig, TopEventInfo& evInfo, edm::ConsumesCollector& iC ):evInfo_(evInfo) {
  typedef std::vector<double> vdouble;

  recoFiltersToken_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));

  csvWeight.initCSVWeight(false, "csvv2");
  bTagWeightL.init(3, "csvv2", BTagEntry::OP_LOOSE , 1);
  bTagWeightM.init(3, "csvv2", BTagEntry::OP_MEDIUM, 1);
  bTagWeightT.init(3, "csvv2", BTagEntry::OP_TIGHT , 1);

  trigTokenMUEL_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = iC.consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));

  jetToken_  = iC.consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = iC.consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = iC.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  const auto muonSet = iConfig.getParameter<edm::ParameterSet>("muon");
  muonToken_ = iC.consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("effSF");
  muonSF_.set(muonSFSet.getParameter<vdouble>("pt_bins"),
              muonSFSet.getParameter<vdouble>("abseta_bins"),
              muonSFSet.getParameter<vdouble>("values"),
              muonSFSet.getParameter<vdouble>("errors"));

  const auto elecSet = iConfig.getParameter<edm::ParameterSet>("electron");
  elecToken_ = iC.consumes<cat::ElectronCollection>(elecSet.getParameter<edm::InputTag>("src"));
  const auto elecSFSet = elecSet.getParameter<edm::ParameterSet>("effSF");
  elecSF_.set(elecSFSet.getParameter<vdouble>("pt_bins"),
              elecSFSet.getParameter<vdouble>("abseta_bins"),
              elecSFSet.getParameter<vdouble>("values"),
              elecSFSet.getParameter<vdouble>("errors"));

  /*
  auto solverPSet = iConfig.getParameter<edm::ParameterSet>("solver");
  auto algoName = solverPSet.getParameter<std::string>("algo");
  std::transform(algoName.begin(), algoName.end(), algoName.begin(), ::toupper);
  if      ( algoName == "CMSKIN" ) solver_.reset(new CMSKinSolver(solverPSet));
  else if ( algoName == "DESYMASSLOOP" ) solver_.reset(new DESYMassLoopSolver(solverPSet));
  else if ( algoName == "DESYSMEARED" ) solver_.reset(new DESYSmearedSolver(solverPSet));
  else if ( algoName == "MT2"    ) solver_.reset(new MT2Solver(solverPSet));
  else if ( algoName == "MAOS"   ) solver_.reset(new MAOSSolver(solverPSet));
  else if ( algoName == "DEFAULT" ) solver_.reset(new TTDileptonSolver(solverPSet));
  else {
    cerr << "The solver name \"" << solverPSet.getParameter<std::string>("algo") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solver_.reset(new TTDileptonSolver(solverPSet)); // A dummy solver
  }
  // PseudoTop
  auto solverPSetPT = iConfig.getParameter<edm::ParameterSet>("solverPseudoTop");
  auto algoNamePT = solverPSetPT.getParameter<std::string>("algo");
  std::transform(algoNamePT.begin(), algoNamePT.end(), algoNamePT.begin(), ::toupper);
  if      ( algoNamePT == "CMSKIN" ) solverPT_.reset(new CMSKinSolver(solverPSetPT));
  else if ( algoNamePT == "DESYMASSLOOP" ) solverPT_.reset(new DESYMassLoopSolver(solverPSetPT));
  else if ( algoNamePT == "DESYSMEARED" ) solverPT_.reset(new DESYSmearedSolver(solverPSetPT));
  else if ( algoNamePT == "MT2"    ) solverPT_.reset(new MT2Solver(solverPSetPT));
  else if ( algoNamePT == "MAOS"   ) solverPT_.reset(new MAOSSolver(solverPSetPT));
  else if ( algoNamePT == "DEFAULT" ) solverPT_.reset(new TTDileptonSolver(solverPSetPT));
  else {
    cerr << "The solver name \"" << solverPSetPT.getParameter<std::string>("algoPT") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solverPT_.reset(new TTDileptonSolver(solverPSetPT)); // A dummy solver
  }
  */
}

void setBranch(TTree* tree,  int sys ) { 

}
 

int TTEventSelector::eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys){
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
  /*
  vector<int> leptonIndex, antiLeptonIndex, jetIndices, bjetIndices;
  VLV allLeptonslv, jetslv;
  vector<double> jetBtags;
  //////////////////////////////////////////////////////// DESY KIN /////////////////////////////////////
  if (selectedBJets.size() > 0){
    LV metlv = evInfo_.metlv;

    int ijet=0;
    for (auto & jet : selectedJets){
      jetslv.push_back(jet.p4());
      jetBtags.push_back(jet.bDiscriminator(BTAG_CSVv2));
      if (jet.bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2L) bjetIndices.push_back(ijet);
      jetIndices.push_back(ijet);
      ++ijet;
    }

    int ilep = 0;
    for (auto & lep : recolep){
      allLeptonslv.push_back(lep->p4());
      if (lep->charge() > 0) antiLeptonIndex.push_back(ilep);
      else leptonIndex.push_back(ilep);
      ++ilep;
    }
    auto kinematicReconstructionSolutions = kinematicReconstruction->solutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices,  allLeptonslv, jetslv, jetBtags, metlv);

    if (evInfo_.step == 5 and sys == sys_nom) evInfo_.cutflow_[10][evInfo_.channel]++;

    if (kinematicReconstructionSolutions.numberOfSolutions()){
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();

      evInfo_.step8 = true;
      if (evInfo_.step == 5)
        if (sys == sys_nom)
          evInfo_.cutflow_[11][evInfo_.channel]++;

      //evInfo_.desytop1 = ToTLorentzVector(top1);
      //evInfo_.desytop2 = ToTLorentzVector(top2);

      LV ttbar = kinematicReconstructionSolutions.solution().ttbar();
      //evInfo_.desyttbar = ToTLorentzVector(ttbar);
      //evInfo_.desyttbar_dphi = deltaPhi(top1.Phi(), top2.Phi());
    }
  }
  */
  return 0;
}

const reco::Candidate* TTEventSelector::getLast(const reco::Candidate* p) const
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

float TTEventSelector::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, sys_e sys) const
{
  float weight = 1.;
  for (auto& m : muons) {
    cat::Muon mu(m);
    if (sys == sys_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == sys_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());

    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.15) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    //weight *= mu.scaleFactor("NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1");
    //weight *= mu.scaleFactor("NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1");
    selmuons.push_back(mu);
  }
  return weight;
}

float TTEventSelector::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
{
  float weight = 1.;
  for (auto& e : elecs) {
    cat::Electron el(e);
    if (sys == sys_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == sys_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) continue;
    //if ( !el.isTrigMVAValid() or !el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") ) continue;
    if (el.relIso(0.3) > 0.12) continue;

    //weight *= el.scaleFactor("mvaEleID-Spring15-25ns-Trig-V1-wp90");
    //weight *= el.scaleFactor("cutBasedElectronID-Spring15-25ns-V1-standalone-medium");
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return weight;
}

cat::JetCollection TTEventSelector::selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, sys_e sys)
{
  // Initialize SF_btag
  float Jet_SF_CSV[19];
  for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] = 1.0;

  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    if (sys == sys_jes_u) jet.setP4(j.p4() * j.shiftedEnUp());
    if (sys == sys_jes_d) jet.setP4(j.p4() * j.shiftedEnDown());
    if (sys == sys_jer_u) jet.setP4(j.p4() * j.smearedResUp());
    if (sys == sys_jer_d) jet.setP4(j.p4() * j.smearedResDown());

    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 2.4)  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep->p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    //if (sys == sys_btag_u) evInfo_.btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 1);
    //else if (sys == sys_btag_d) evInfo_.btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, -1);
    //else evInfo_.btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 0);
    for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] *= csvWeight.getSF(jet, iu);
    seljets.push_back(jet);
  }
  for (unsigned int iu=0; iu<19; iu++) evInfo_.csvweights.push_back(Jet_SF_CSV[iu]);

  evInfo_.btagweight = Jet_SF_CSV[0];
  // if      ( sys == sys_btag_u ) evInfo_.btagweight = bTagWeightL.eventWeight(seljets, 1);
  // else if ( sys == sys_btag_d ) evInfo_.btagweight = bTagWeightL.eventWeight(seljets, 2);
  // else                          evInfo_.btagweight = bTagWeightL.eventWeight(seljets, 0);

  return seljets;
}

cat::JetCollection TTEventSelector::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2L) continue;
    //if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;//forsync
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}
 
