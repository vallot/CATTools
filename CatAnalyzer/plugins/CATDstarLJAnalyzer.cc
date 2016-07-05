#include "CATTools/CatAnalyzer/interface/CATDstarAnalyzer.h"

using namespace std;
using namespace cat;
using namespace dileptonCommonGlobal;

class CATDstarLJAnalyzer : virtual public CATDstarAnalyzer {
  public:
    explicit CATDstarLJAnalyzer(const edm::ParameterSet& iConfig);
    ~CATDstarLJAnalyzer();
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) override;
    virtual float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, cat::MuonCollection& vetomuons, sys_e sys) const;
    virtual float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, cat::ElectronCollection& vetoelecs, sys_e sys) const;
    bool isSelectMuon(cat::Muon) const ;
    bool isVetoMuon  (cat::Muon) const ;
    bool isSelectElec(cat::Electron) const ;
    bool isVetoElec  (cat::Electron) const ;

    virtual cat::JetCollection selectBJets(const JetCollection& jets) const override;
    
  protected:
    edm::EDGetTokenT<int> trigTokenMUJET_;
    edm::EDGetTokenT<int> trigTokenELJET_;

  private:
};

CATDstarLJAnalyzer::~CATDstarLJAnalyzer(){
  cout <<setw(10)<<"cut flow"<<setw(10)<<"no Lep+Jet"<<setw(10)<<"mu+jet"<<setw(10)<<"e+jet"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<setw(10)<<"step "<<i<< setw(10)<<cutflow_[i][0] << setw(10)<<cutflow_[i][1] << setw(10)<<cutflow_[i][2] << endl;
  }
}

bool CATDstarLJAnalyzer::isVetoMuon(cat::Muon mu) const
{
  if (!mu.isLooseMuon() ) return false;
  if (mu.pt() <10. ) return false;
  if (std::abs(mu.eta())>2.4 ) return false;
  if (mu.relIso(0.4) > 0.25 ) return false;
  return true; 
}
bool CATDstarLJAnalyzer::isSelectMuon(cat::Muon mu) const
{
  if (!mu.isTightMuon() ) return false;
  if (mu.pt() <26. ) return false;
  if (std::abs(mu.eta())>2.1 ) return false;
  if (mu.relIso(0.4) > 0.15 ) return false;
  return true; 
}

float CATDstarLJAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, cat::MuonCollection& vetomuons, sys_e sys ) const
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

bool CATDstarLJAnalyzer::isVetoElec(cat::Electron el) const
{
  if (el.pt() < 15. ) return false;
  if (std::abs(el.eta()) > 2.4) return false;
  if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) return false;

  /*
  if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) return false;
  if ( std::abs(el.scEta()) <= 1.479 ) { if ( el.relIso(0.3) > 0.126) return false; }
  else                                 { if ( el.relIso(0.3) > 0.144) return false; }
  */

  return true;
}

bool CATDstarLJAnalyzer::isSelectElec(cat::Electron el) const
{
  if (el.pt() < 30. ) return false;
  if (std::abs(el.eta()) > 2.4) return false;
  if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) return false;

  /*
  if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) return false;
  if ( std::abs(el.scEta()) <= 1.479 ) { if ( el.relIso(0.3) > 0.0766) return false; }
  else                                 { if ( el.relIso(0.3) > 0.0678) return false; }
  */

  return true;
}

float CATDstarLJAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, cat::ElectronCollection& vetoelecs, sys_e sys) const
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

CATDstarLJAnalyzer::CATDstarLJAnalyzer(const edm::ParameterSet& iConfig) : CATDstarAnalyzer(iConfig)
{
  trigTokenMUJET_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUJET"));
  trigTokenELJET_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELJET"));
}

void CATDstarLJAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  const bool runOnMC = !iEvent.isRealData();
  //cutflow_[0][0]++;

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();
    if( sys == 0 ) genInfo(iEvent, iSetup);
    int terminate = eventSelection(iEvent, iSetup,sys);
    if ( terminate == -1 ) continue;
    else if ( terminate == -2 ) return;

    analyzeCustom(iEvent, iSetup, sys);
    ttree_[sys]->Fill();
  }
}

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
  if (sys == sys_nom) cutflow_[0][b_channel]++;

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
  cat::MuonCollection selMuons , vetoMuons;
  cat::ElectronCollection selElecs, vetoElecs;
  selectMuons(*muons, selMuons, vetoMuons, (sys_e)sys);
  selectElecs(*electrons, selElecs, vetoElecs, (sys_e)sys);
  if ( selMuons.size()+selElecs.size() ==0 ) {
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

  const cat::Lepton& recolep1 = *recolep[0];


  // Determine channel
  const int pdgIdSum = std::abs(recolep1.pdgId()) ;
  if (pdgIdSum == 13) b_channel = CH_MUJET; // mu+jet
  if (pdgIdSum == 11) b_channel = CH_ELJET; // el+jet
  /*
  if      ( selElecs.size() > 0)  b_channel = CH_ELJET;
  else if ( selMuons.size() > 0 ) b_channel = CH_MUJET;
  */

  b_mueffweight    = getMuEffSF(recolep1,  0);
  b_mueffweight_up = getMuEffSF(recolep1, +1);
  b_mueffweight_dn = getMuEffSF(recolep1, -1);

  b_eleffweight    = getElEffSF(recolep1,  0);
  b_eleffweight_up = getElEffSF(recolep1, +1);
  b_eleffweight_dn = getElEffSF(recolep1, -1);

  // Trigger results
  // No information for Single lepton trigger
  b_tri = 0;

  edm::Handle<int> trigHandle;
  if      ( b_channel == CH_ELJET ) iEvent.getByToken(trigTokenELJET_, trigHandle);
  else if ( b_channel == CH_MUJET ) iEvent.getByToken(trigTokenMUJET_, trigHandle);
  if ( *trigHandle != 0 ){
    //std::cout<<"Active trigger!"<<std::endl;
    b_tri = 1.0; // Need to Update
    b_tri_up = 1.0;
    b_tri_dn = 1.0;
  }
  //else  std::cout<<"No trigger!"<<std::endl;

  b_lep1 = recolep1.tlv(); b_lep1_pid = recolep1.pdgId();
  //b_lep2 = recolep2.tlv(); b_lep2_pid = recolep2.pdgId();
  b_lep2 = TLorentzVector(); b_lep2_pid = 0;
  b_dilep = TLorentzVector();

  b_step1 = true;
  b_step = 1;
  if (sys == sys_nom) cutflow_[4][b_channel]++;

  if ( b_channel == CH_MUJET && vetoMuons.size() ==0 ) {
    if ( b_step1 ) {
      b_step2 = true;
      b_step=2 ;
      if (sys == sys_nom) cutflow_[5][b_channel]++;
    }
  }
  if ( b_channel == CH_ELJET && vetoElecs.size() ==0 ) {
    if ( b_step1 ) {
      b_step2 = true;
      b_step=2 ;
      if (sys == sys_nom) cutflow_[5][b_channel]++;
    }
  }

  selectedJets.clear();
  selectedBJets.clear();

  selectedJets  = selectJets(*jets, recolep, (sys_e)sys);
  selectedBJets =  selectBJets(selectedJets);

  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();

  if ( b_njet >=1 ){
    if (b_step2){
      b_step3 = true;
      ++b_step;
      if (sys == sys_nom) cutflow_[6][b_channel]++;
    }
  }
  if ( b_njet >=4 ){
    if (b_step3){
      b_step4 = true;
      ++b_step;
      if (sys == sys_nom) cutflow_[7][b_channel]++;
    }
  }
  if ( b_nbjet >= 1){
    if (b_step4){
      b_step5 = true;
      ++b_step;
      if (sys == sys_nom) cutflow_[8][b_channel]++;
    }
  }
  if ( b_nbjet >= 2){
    if (b_step5){
      b_step6 = true;
      ++b_step;
      if (sys == sys_nom) cutflow_[9][b_channel]++;
    }
  }
  if ( b_nbjet == 2){
    if (b_step6){
      b_step7 = true;
      ++b_step;
      if (sys == sys_nom) cutflow_[10][b_channel]++;
    }
  }

  return 0;
}

cat::JetCollection CATDstarLJAnalyzer::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;//forsync
      selBjets.push_back(jet);
    }
    return selBjets;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATDstarLJAnalyzer);
