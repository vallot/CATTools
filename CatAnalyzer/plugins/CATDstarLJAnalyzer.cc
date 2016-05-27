#include"CATDstarAnalyzer.h"


class CATDstarLJAnalyzer : public CATDstarAnalyzer {
  public:
    explicit CATDstarLJAnalyzer(const edm::ParameterSet&);
    ~CATDstarLJAnalyzer();

  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) override;
};

int CATDstarLJAnalyzer::eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys){
  const bool runOnMC = !iEvent.isRealData();

  if (runOnMC){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    b_puweight = *puweightHandle;

    edm::Handle<float> puweightHandle_up;
    iEvent.getByToken(puweightToken_up_, puweightHandle_up);
    b_puweight_up = *puweightHandle_up;

    edm::Handle<float> puweightHandle_dn;
    iEvent.getByToken(puweightToken_dn_, puweightHandle_dn);
    b_puweight_dn = *puweightHandle_dn;

    edm::Handle<float> genweightHandle;
    iEvent.getByToken(genWeightToken_, genweightHandle);
    b_genweight = (*genweightHandle);
    b_weight = b_genweight*b_puweight;

  }

  if (sys == sys_nom) h_nevents->Fill(0.5,b_puweight*b_genweight);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    if (b_keepTtbarSignal) ttree_[sys]->Fill();
    std::cout<<"No PV"<<std::endl;
    return -1;
  }
  if (sys == sys_nom) cutflow_[1][b_channel]++;

  // const reco::Vertex &PV = vertices->front();
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  b_nvertex = *nGoodVertexHandle;

  edm::Handle<int> lumiSelectionHandle;
  iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
  if (!runOnMC){
    if (*lumiSelectionHandle == 0) return -2;  // Critical problem! Terminate the job for this sample.
  }

  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;
  if (sys == sys_nom) cutflow_[2][b_channel]++;

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
    if (b_keepTtbarSignal) ttree_[sys]->Fill();
    return -1;
  }
  if (sys == sys_nom) cutflow_[3][b_channel]++;

  //std::vector<const cat::Lepton*> recolep;
  auto& recolep = recolep_;
  recolep.clear();
  for ( const auto& x : selMuons ) recolep.push_back(&x);
  for ( const auto& x : selElecs ) recolep.push_back(&x);

  sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
  b_is3lep = recolep.size();
  recolep.erase(recolep.begin()+2,recolep.end());
  if ( recolep[0]->charge() < 0 ){
    swap(recolep[0], recolep[1]);
  }
  const cat::Lepton& recolep1 = *recolep[0];
  const cat::Lepton& recolep2 = *recolep[1];

  // Determine channel
  const int pdgIdSum = std::abs(recolep1.pdgId()) + std::abs(recolep2.pdgId());
  if (pdgIdSum == 24) b_channel = CH_MUEL; // emu
  if (pdgIdSum == 22) b_channel = CH_ELEL; // ee
  if (pdgIdSum == 26) b_channel = CH_MUMU; // mumu

  b_mueffweight    = getMuEffSF(recolep1,  0)*getMuEffSF(recolep2,  0);
  b_mueffweight_up = getMuEffSF(recolep1, +1)*getMuEffSF(recolep2, +1);
  b_mueffweight_dn = getMuEffSF(recolep1, -1)*getMuEffSF(recolep2, -1);

  b_eleffweight    = getElEffSF(recolep1,  0)*getElEffSF(recolep2,  0);
  b_eleffweight_up = getElEffSF(recolep1, +1)*getElEffSF(recolep2, +1);
  b_eleffweight_dn = getElEffSF(recolep1, -1)*getElEffSF(recolep2, -1);

  // Trigger results
  // Scale factors are from AN16-025 (v4) http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_025_v4.pdf
  b_tri = 0;
  edm::Handle<int> trigHandle;
  if ( b_channel == CH_ELEL ) {
    iEvent.getByToken(trigTokenELEL_, trigHandle);
    if ( *trigHandle != 0 ) b_tri = 0.953; // +- 0.009
  }
  else if ( b_channel == CH_MUMU ) {
    iEvent.getByToken(trigTokenMUMU_, trigHandle);
    if ( *trigHandle != 0 ) b_tri = 0.948; // +- 0.002
  }
  else if ( b_channel == CH_MUEL ) {
    iEvent.getByToken(trigTokenMUEL_, trigHandle);
    if ( *trigHandle != 0 ) b_tri = 0.975; // +- 0.004
  }

  b_lep1 = recolep1.tlv(); b_lep1_pid = recolep1.pdgId();
  b_lep2 = recolep2.tlv(); b_lep2_pid = recolep2.pdgId();
  b_dilep = b_lep1+b_lep2;
  const auto tlv_ll = recolep1.p4()+recolep2.p4();

  if (tlv_ll.M() < 20. || recolep1.charge() * recolep2.charge() > 0){
    if (b_keepTtbarSignal) ttree_[sys]->Fill();
    return -1;
  }
  b_step1 = true;
  b_step = 1;
  if (sys == sys_nom) cutflow_[4][b_channel]++;

  if ( (b_channel == CH_MUEL) || ((tlv_ll.M() < 76) || (tlv_ll.M() > 106)) ){
    b_step2 = true;
    b_step = 2;
    if (sys == sys_nom) cutflow_[5][b_channel]++;
  }

  selectedJets.clear();
  selectedBJets.clear();

  selectedJets  = selectJets(*jets, recolep, (sys_e)sys);
  selectedBJets =  selectBJets(selectedJets);

  //met.clear();
  met = mets->front().p4();
  b_met = met.pt();
  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();

  if (selectedJets.size() >1 ){
    b_step3 = true;
    if (b_step == 2){
      ++b_step;
      if (sys == sys_nom) cutflow_[6][b_channel]++;
    }
  }

  if ((b_channel == CH_MUEL) || (b_met > 40.)){
    b_step4 = true;
    if (b_step == 3){
      ++b_step;
      if (sys == sys_nom) cutflow_[7][b_channel]++;
    }
  }

  if (selectedBJets.size() > 0){
    b_step5 = true;
    if (b_step == 4){
      ++b_step;
      if (sys == sys_nom) cutflow_[8][b_channel]++;
    }
  }
  vector<int> leptonIndex, antiLeptonIndex, jetIndices, bjetIndices;
  VLV allLeptonslv, jetslv;
  vector<double> jetBtags;
  //////////////////////////////////////////////////////// DESY KIN /////////////////////////////////////
  if (selectedBJets.size() > 0){
    LV metlv = mets->front().p4();

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
    KinematicReconstructionSolutions kinematicReconstructionSolutions  =  kinematicReconstruction->solutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices,  allLeptonslv, jetslv, jetBtags, metlv);

    if (b_step == 5 and sys == sys_nom) cutflow_[10][b_channel]++;

    if (kinematicReconstructionSolutions.numberOfSolutions()){
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();

      b_step8 = true;
      if (b_step == 5)
        if (sys == sys_nom)
          cutflow_[11][b_channel]++;

      b_desytop1 = ToTLorentzVector(top1);
      b_desytop2 = ToTLorentzVector(top2);

      LV ttbar = kinematicReconstructionSolutions.solution().ttbar();
      b_desyttbar = ToTLorentzVector(ttbar);
      b_desyttbar_dphi = deltaPhi(top1.Phi(), top2.Phi());
    }
  }
  return 0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(CATDstarLJAnalyzer);
