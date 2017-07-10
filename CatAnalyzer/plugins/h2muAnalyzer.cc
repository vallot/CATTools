/////////////////////////////////////////////////////////////////////////////////////////////////////
//only jets that pass though the jet selection are considered for the b-jet veto 

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenWeights.h"

#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "CATTools/CatAnalyzer/interface/analysisUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TTree.h"
#include "TH1D.h"
#include "CATTools/CatAnalyzer/src/RoccoR.cc"
//#include "CATTools/CatAnalyzer/src/RoccoR.h"
//#include "TLorentzVector.h"

using namespace std;
using namespace cat;

class h2muAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit h2muAnalyzer(const edm::ParameterSet&);
  ~h2muAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override{};
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  void setBranch(TTree* tree, systematic sys);
  void resetBranch();
  bool eventSelection(const edm::Event& iEvent, systematic sys);
  
  MuonCollection selectMuons(const MuonCollection& muons, systematic sys);
  ElectronCollection selectElectrons(const ElectronCollection& elecs, MuonCollection& recolep, systematic sys);
  JetCollection selectJets(const JetCollection& jets, MuonCollection& recolep, systematic sys);
  JetCollection selectBJets(const JetCollection& jets);
  bool preSelection(const JetCollection& seljets, float MET) const;
  int jetCategory(const JetCollection& seljets, float MET, float ll_pt) const;
  int etaCategory(float lep1_eta, float lep2_eta) const;
  bool bumpCategory(const JetCollection& jets, const MuonCollection& recolep);
  bool bumpCat1(const JetCollection& jets, const MuonCollection& recolep);
  bool bumpCat2(const JetCollection& jets, const MuonCollection& recolep, float met);
  ScaleFactorEvaluator muonSF_, elecSF_;
  
  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genweightToken_, puweightToken_, puweightToken_up_, puweightToken_dn_, topPtWeight_;
  //edm::EDGetTokenT<vector<float>> pdfweightToken_, scaleweightToken_;
  edm::EDGetTokenT<MuonCollection>     muonToken_;
  edm::EDGetTokenT<ElectronCollection> elecToken_;
  edm::EDGetTokenT<JetCollection>      jetToken_;
  edm::EDGetTokenT<METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;  
  //edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  vector<edm::EDGetTokenT<edm::TriggerResults>> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  RoccoR *rocCor;

  bool runOnMC;
  vector<TTree*> ttree_;
  TH1D *h_nevents, *h_cutflow, *h_cutflowCat, *h_Tgenweight ;
  int b_run, b_lumi, b_event;
  int b_nvertex, b_step, b_channel, b_njet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_filtered;
  float b_tri;
  float b_met, b_weight, b_puweight, b_puweight_up, b_puweight_dn, b_genweight,
    b_mueffweight, b_mueffweight_up, b_mueffweight_dn,b_beffweight, b_trigeffweight,
    b_eleffweight, b_eleffweight_up, b_eleffweight_dn, b_genweightT;
  vector<float> b_pdfWeights, b_scaleWeights;

  int size;
  int b_cat;
  int b_cat_eta;
  bool b_isLoose, b_isMedium, b_isTight;
  
  double scaleFactor;
  double u1, u2;
  bool b_bumpcat1, b_bumpcat2;
  int event;
  TLorentzVector b_genlep1; int b_genlep1_pid;
  TLorentzVector b_genlep2; int b_genlep2_pid;
  TLorentzVector b_gendilep;
  
  TLorentzVector b_lep1; int b_lep1_pid;
  TLorentzVector b_lep2; int b_lep2_pid;
  TLorentzVector b_jet1, b_jet2, b_dilep, b_dijet;

};

h2muAnalyzer::h2muAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  //pdfweightToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("pdfweight"));
  //scaleweightToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaleweight"));
  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightToken_up_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_up"));
  puweightToken_dn_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_dn"));
  jetToken_  = consumes<JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));
  //triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  for ( auto x : iConfig.getParameter<vector<edm::InputTag>>("triggerBits") ) {
    triggerBits_.push_back(consumes<edm::TriggerResults>(x));
  }
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  typedef vector<double> vdouble;

  const auto muonSet = iConfig.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("effSF");
  muonSF_.set(muonSFSet.getParameter<vdouble>("pt_bins"),
              muonSFSet.getParameter<vdouble>("eta_bins"),
              muonSFSet.getParameter<vdouble>("values"),
              muonSFSet.getParameter<vdouble>("errors"));

  const auto elecSet = iConfig.getParameter<edm::ParameterSet>("electron");
  elecToken_ = consumes<ElectronCollection>(elecSet.getParameter<edm::InputTag>("src"));
  const auto elecSFSet = elecSet.getParameter<edm::ParameterSet>("effSF");
  elecSF_.set(elecSFSet.getParameter<vdouble>("pt_bins"),
              elecSFSet.getParameter<vdouble>("eta_bins"),
              elecSFSet.getParameter<vdouble>("values"),
              elecSFSet.getParameter<vdouble>("errors"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_nevents = fs->make<TH1D>("nevents","nevents",4,0,4);       
  h_cutflow = fs->make<TH1D>("cutflow","cutflow",15,0,15);       
  h_cutflowCat = fs->make<TH1D>("cutflowCat","cutflowCat",6,0,6);     
  for (int sys = 0; sys < syst_total; ++sys){
    ttree_.push_back(fs->make<TTree>(systematicName[systematic(sys)].c_str(), systematicName[systematic(sys)].c_str()));
    auto tr = ttree_.back();
    setBranch(tr, systematic(sys));
  }
  
 //rocCor = new RoccoR(edm::FileInPath("CATTools/CatAnalyzer/data/rcdata.2016.v3/").fullPath());  
 rocCor = new RoccoR(std::string(std::getenv("CMSSW_BASE"))+"/src/CATTools/CatAnalyzer/data/rcdata.2016.v3/");

}

h2muAnalyzer::~h2muAnalyzer()
{
  cout <<"     cut flow "<< endl;
  for ( int i=1; i<12; ++i ) {
    cout <<"step "<< i << "    "<< h_cutflow->GetBinContent(i+1) << endl;
  }
  for ( int i=1; i<6; ++i ) {
    cout <<"Categories "<< i << "    "<< h_cutflowCat->GetBinContent(i+1) << endl;
  }  
}

void h2muAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC = !iEvent.isRealData();
  bool keepEvent = false;
  
  for (int sys = 0; sys < syst_total; ++sys){
    if (sys != syst_nom && !runOnMC) break;
    resetBranch();
    
    keepEvent = eventSelection(iEvent, systematic(sys));
    
    if (keepEvent)
      ttree_[sys]->Fill();
    break;
  }
}

bool h2muAnalyzer::eventSelection(const edm::Event& iEvent, systematic sys)
{
  if (sys == syst_nom) h_cutflow->Fill(0);
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();
  if (runOnMC && sys == syst_nom){
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
    iEvent.getByToken(genweightToken_, genweightHandle);
    b_genweight = (*genweightHandle);
    b_weight = b_genweight*b_puweight;
    b_trigeffweight = 1;
    b_beffweight = 1;
    h_nevents->Fill(0.5,b_weight);
    h_nevents->Fill(2.5,b_genweight);
    h_nevents->Fill(3.5,1);

    edm::Handle<reco::GenParticleCollection> genParticles;
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
      b_genlep1 = genLep1;
      b_genlep2 = genLep2;
      b_genlep1_pid = -13;
      b_genlep2_pid = 13;
      b_gendilep = b_genlep1 + b_genlep2;
    }
  }
  //if (b_event != 112520) return false; 
  edm::Handle<int> lumiSelectionHandle;
  iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
  if (!runOnMC){
    if (*lumiSelectionHandle == 0) return false;
  }
  b_step1 = true;
  b_step = 1;
  if (sys == syst_nom) h_cutflow->Fill(1);
 
  ////////////////////////////////////////////////////////////////////////////////
  // filters
  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;

  ////////////////////////////////////////////////////////////////////////////////
  // vertex
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    return false;
  }  
  b_step2 = true;
  b_step = 2;
  if (sys == syst_nom) h_cutflow->Fill(2);

  // const reco::Vertex &PV = vertices->front();
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  b_nvertex = *nGoodVertexHandle;

  ////////////////////////////////////////////////////////////////////////////////
  // muon selection
  edm::Handle<MuonCollection> muons; iEvent.getByToken(muonToken_, muons);
  MuonCollection selectedMuons = selectMuons(*muons, sys);
  if ( selectedMuons.size() != 2 ) return false;
  if (sys == syst_nom) h_cutflow->Fill(3);

  const Muon &mu1 = selectedMuons[0], &mu2 = selectedMuons[1];
  if (mu1.charge()*mu2.charge() > 0) return false;
  if (sys == syst_nom) h_cutflow->Fill(4);
  b_step3 = true;
  b_step = 3;

  b_lep1 = mu1.tlv();b_lep1_pid = mu1.pdgId();
  b_lep2 = mu2.tlv();b_lep2_pid = mu2.pdgId();
  b_dilep = b_lep1 + b_lep2;

  if (b_dilep.M() < 12) return false; 
  if (sys == syst_nom) h_cutflow->Fill(4);
  b_isLoose = (mu1.isLooseMuon() && mu2.isLooseMuon());
  b_isMedium = (mu1.isMediumMuon() && mu2.isMediumMuon());
  b_isTight = (mu1.isTightMuon() && mu2.isTightMuon());
    
   b_mueffweight    = muonSF_.getScaleFactor(mu1, 13, 0)*muonSF_.getScaleFactor(mu2, 13,  0);
   b_mueffweight_up = muonSF_.getScaleFactor(mu1, 13, +1)*muonSF_.getScaleFactor(mu2, 13, +1);
   b_mueffweight_dn = muonSF_.getScaleFactor(mu1, 13, -1)*muonSF_.getScaleFactor(mu2, 13, -1);
  
  ///////////////////////////////////////////////////////////////////////////////
  //electron veto
  edm::Handle<ElectronCollection> elecs; iEvent.getByToken(elecToken_, elecs);
  ElectronCollection selectedElecs = selectElectrons(*elecs, selectedMuons, sys);
  if (selectedElecs.size() > 0) return true;
  if (sys == syst_nom) h_cutflow->Fill(6);
  
  ////////////////////////////////////////////////////////////////////////////////
  //b-jet veto
  edm::Handle<JetCollection> jets; iEvent.getByToken(jetToken_, jets);
  JetCollection selectedJets = selectJets(*jets, selectedMuons, sys);
  JetCollection selectedBJets = selectBJets(selectedJets);
      
  if (selectedBJets.size() > 0 ) return true;
  if (sys == syst_nom) h_cutflow->Fill(7);
  
  ////////////////////////////////////////////////////////////////////////////////
  // trigger
  edm::Handle<edm::TriggerResults> triggerBits;
  for ( auto token : triggerBits_ ) {
    if ( iEvent.getByToken(token, triggerBits) ) break;
  }
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  if (!trigHelper.triggerFired("HLT_IsoMu24_v") && !trigHelper.triggerFired("HLT_IsoTkMu24_v")){
    return true;  
  }
  if (sys == syst_nom) h_cutflow->Fill(8);
  b_step4 = true;
  b_step = 4;  
  
  ////////////////////////////////////////////////////////////////////////////////
  // trigger matching
  bool triggerMatching = false;
  if ( mu1.pt() > 26 &&
       (trigHelper.triggerMatched("HLT_IsoMu24_v", mu1) ||
        trigHelper.triggerMatched("HLT_IsoTkMu24_v", mu1)))
    triggerMatching = true;
    
  if ( mu2.pt() > 26 &&
       (trigHelper.triggerMatched("HLT_IsoMu24_v", mu2) ||
        trigHelper.triggerMatched("HLT_IsoTkMu24_v", mu2)))
    triggerMatching = true;

  if (!triggerMatching) return true;                            
  if (sys == syst_nom) h_cutflow->Fill(9);
  b_step5 = true;
  b_step = 5;

  
  ////////////////////////////////////////////////////////////////////////////////
  // jets  
  b_njet = selectedJets.size();
  if (b_njet > 0)
    b_jet1 = selectedJets[0].tlv();
  if (b_njet > 1){
    b_jet2 = selectedJets[1].tlv();
    b_dijet = b_jet1 + b_jet2;
  }

  edm::Handle<METCollection> mets; iEvent.getByToken(metToken_, mets);
  const auto met = mets->front().p4();
  b_met = met.pt();
  b_cat_eta = etaCategory(b_lep1.Eta(), b_lep2.Eta());
  b_cat = jetCategory(selectedJets, b_met, b_dilep.Pt());
  h_cutflowCat->Fill(b_cat);
  h_nevents->Fill(1.5,b_genweight);
  ofstream eventdump;
  eventdump.open("Event_Dump.txt", ios::app);
  for (int i = 0; i < b_njet; i++)
    {
    eventdump<< "check jets1 "<<b_event<< " jet"<<i<<"->pt():" << selectedJets[i].tlv().Pt() << " jet"<<i<<"->eta():" << selectedJets[i].tlv().Eta() <<" jet"<<i<<"->phi():" << selectedJets[i].tlv().Phi()<<endl;
    }
    eventdump<<""<<endl;
  eventdump.close();

  return true;
}
MuonCollection h2muAnalyzer::selectMuons(const MuonCollection& muons, systematic sys)
{
  MuonCollection selmuons;
  for (auto& m : muons) {
    Muon mu(m);
    if (!mu.isGlobalMuon()) continue;
    if (!mu.isTrackerMuon()) continue;
    if (abs(mu.eta()) > 2.4) continue;
    if (!mu.isMediumMuon()) continue;  //<------------!

    // TLorentzVector tmu(mu.tlv());
    // /*
    //   cout << "\n initial   " << mu.tlv().Pt()
    //   << ", " << mu.tlv().Eta()
    //   << ", " << mu.tlv().Phi()
    //   << ", " << mu.tlv().M()
    //   << " mu.trackerLayersWithMeasurement() " << mu.trackerLayersWithMeasurement()
    //   << " mu.numberOfValidHits() " << mu.numberOfValidHits()
    //   <<endl;
    // */
    scaleFactor = 1.0;
    u1 = gRandom->Rndm();
    u2 = gRandom->Rndm();
    if (!runOnMC){ 
        scaleFactor = rocCor -> kScaleDT(mu.charge(), mu.pt(), mu.eta(), mu.phi(), 0, 0);
    }
    else{
        if (mu.pt() == b_genlep1.Pt())
            scaleFactor = rocCor -> kScaleFromGenMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), mu.trackerLayersWithMeasurement(), b_genlep1.Pt(), u1, 0, 0);
        if (mu.pt() == b_genlep2.Pt())
            scaleFactor = rocCor -> kScaleFromGenMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), mu.trackerLayersWithMeasurement(), b_genlep2.Pt(), u1, 0, 0);
        else {
            scaleFactor = rocCor -> kScaleAndSmearMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), mu.trackerLayersWithMeasurement(), u1, u2, 0, 0);
        }
    } 
    mu.setP4(m.p4() * scaleFactor);

    //cout << " before   " << mu.tlv().Pt()
    //<< ", " << mu.tlv().Eta()
    //<< ", " << mu.tlv().Phi()
    //<< ", " << mu.tlv().M()
    //<<endl;

    //cout << "initial pt:     "<<mu.tlv().Pt()<< endl;
    
    //cout << " applied1   " << mu.pt()
    //<< ", " << mu.tlv().Eta()
    //<< ", " << mu.tlv().Phi()
    //<< ", " << mu.tlv().M()
    //<<endl;

    // float qter = 1.0;
    // if (runOnMC)
    //   rocCor->momcor_mc(tmu, mu.charge(), mu.trackerLayersWithMeasurement(), qter);
    // else 
    //   rocCor->momcor_data(tmu, mu.charge(), 0, qter);

    // mu.setP4(m.p4() * tmu.E()/mu.tlv().E());
    /*
      if (tmu.Pt() != mu.tlv().Pt()){
      cout << " corrected " << tmu.Pt()
      << ", " << tmu.Eta()
      << ", " << tmu.Phi()
      << ", " << tmu.M()
      <<endl;
    
      cout << " applied   " << mu.tlv().Pt()
      << ", " << mu.tlv().Eta()
      << ", " << mu.tlv().Phi()
      << ", " << mu.tlv().M()
      <<endl;
      }
    */
    if (sys == syst_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == syst_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());
    if (mu.pt() < 10.) continue;
    if (mu.pt() > 20.) continue;
    if (mu.relIso(0.4) > 0.25) continue;
    selmuons.push_back(mu);
  }
  return selmuons;
}

ElectronCollection h2muAnalyzer::selectElectrons(const ElectronCollection& elecs, MuonCollection& recolep, systematic sys)
{
  ElectronCollection selelecs;
  for (auto& e : elecs) {
    Electron el(e);
    if (abs(el.scEta()) > 2.5) continue;
    if ((abs(el.scEta()) > 1.4442) && (abs(el.scEta()) < 1.566)) continue;
    if (!el.electronID("cutBasedElectronID-Summer16-80X-V1-medium") ) continue;
    if (sys == syst_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == syst_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if (el.pt() < 10.) continue;
    if (el.relIso(0.3) > 0.15) continue;
    bool hasOverLap1 = false;
    for (auto muon : recolep){
      if (deltaR(el.p4(),muon.p4()) < 0.4) hasOverLap1 = true;
    }
    if (hasOverLap1) continue;
    selelecs.push_back(el);
  }
  return selelecs;
}

JetCollection h2muAnalyzer::selectJets(const JetCollection& jets, MuonCollection& recolep, systematic sys)
{
  JetCollection seljets;
  for (auto& j : jets) {
    Jet jet(j);
    if (abs(jet.eta()) > 4.7) continue;
    if (!jet.LooseId()) continue;
    
    if (sys == syst_jes_u) jet.setP4(j.p4() * j.shiftedEnUp());
    if (sys == syst_jes_d) jet.setP4(j.p4() * j.shiftedEnDown());
    if (sys == syst_jer_u) jet.setP4(j.p4() * j.smearedResUp());
    if (sys == syst_jer_d) jet.setP4(j.p4() * j.smearedResDown());

    if (jet.pt() < 30.) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep.p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    seljets.push_back(jet);
  }
  return seljets;
}

JetCollection h2muAnalyzer::selectBJets(const JetCollection& jets)
{
    JetCollection selBjets;
    for (auto& jet : jets) { 
        //if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2L) continue;
        if (abs(jet.eta()) > 2.4)  continue; 
        if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;//forsync
        if (jet.pt() < 30.) continue;                                    
        //printf("b jet with pt %4.1f\n", jet.pt());
        selBjets.push_back(jet);
    }
    return selBjets;
}
void h2muAnalyzer::setBranch(TTree* tr, systematic sys)
{
  tr->Branch("step", &b_step, "step/I");
  tr->Branch("channel", &b_channel, "channel/I");
  tr->Branch("step1", &b_step1, "step1/O");
  tr->Branch("step2", &b_step2, "step2/O");
  tr->Branch("step3", &b_step3, "step3/O");
  tr->Branch("step4", &b_step4, "step4/O");
  tr->Branch("step5", &b_step5, "step5/O");
  tr->Branch("step6", &b_step5, "step6/O");
  tr->Branch("tri", &b_tri, "tri/F");
  tr->Branch("filtered", &b_filtered, "filtered/O");
    
  tr->Branch("nvertex", &b_nvertex, "nvertex/I");
  tr->Branch("njet", &b_njet, "njet/I");
  tr->Branch("met", &b_met, "met/F");
  tr->Branch("isLoose", &b_isLoose, "isLoose/O");
  tr->Branch("isMedium", &b_isMedium, "isMedium/O");
  tr->Branch("isTight", &b_isTight, "isTight/O");
  tr->Branch("lep1", "TLorentzVector", &b_lep1);
  tr->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
  tr->Branch("lep2", "TLorentzVector", &b_lep2);
  tr->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
  tr->Branch("dilep", "TLorentzVector", &b_dilep);
  tr->Branch("jet1", "TLorentzVector", &b_jet1);
  tr->Branch("jet2", "TLorentzVector", &b_jet2);
  tr->Branch("dijet", "TLorentzVector", &b_dijet);
  
  //final hierachy
  //(e.g. In case of 0,1jet, Tight and Loose.Otherwise 2jet include VBF Tight, ggF Tight, Loose)
  tr->Branch("cat", &b_cat, "cat/I");
  //Geometrical Categorization
  //only included 0jet and 1jet
  tr->Branch("cat_eta", &b_cat_eta, "cat_eta/I");
  
  //28GeV bump category
  tr->Branch("bumpcat1", &b_bumpcat1, "bumpcat1/O");
  tr->Branch("bumpcat2", &b_bumpcat2, "bumpcat2/O");

  tr->Branch("weight", &b_weight, "weight/F");
  tr->Branch("puweight", &b_puweight, "puweight/F");
  tr->Branch("genweight", &b_genweight, "genweight/F");
  tr->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  tr->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  tr->Branch("trigeffweight", &b_trigeffweight, "trigeffweight/F");
  tr->Branch("beffweight", &b_beffweight, "beffweight/F");
  tr->Branch("genweightT",&b_genweightT, "genweightT/F");
  // only save for nomial ttree 
  if (sys != syst_nom) return;
  tr->Branch("run", &b_run, "run/I");
  tr->Branch("event", &b_event, "event/I");
    
  tr->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
  tr->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
  tr->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  tr->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
  tr->Branch("eleffweight_up", &b_eleffweight_up, "eleffweight_up/F");
  tr->Branch("eleffweight_dn", &b_eleffweight_dn, "eleffweight_dn/F");      
  tr->Branch("pdfWeights","vector<float>",&b_pdfWeights);
  tr->Branch("scaleWeights","vector<float>",&b_scaleWeights);
      
  tr->Branch("genlep1", "TLorentzVector", &b_genlep1);
  tr->Branch("genlep1_pid", &b_genlep1_pid, "genlep1_pid/I");    
  tr->Branch("genlep2", "TLorentzVector", &b_genlep2);
  tr->Branch("genlep2_pid", &b_genlep2_pid, "genlep2_pid/I");    
  tr->Branch("gendilep", "TLorentzVector", &b_gendilep);  
}

void h2muAnalyzer::resetBranch()
{
  b_nvertex = 0;b_step = 0;b_channel = 0;b_njet = 0;
  b_step1 = 0;b_step2 = 0;b_step3 = 0;b_step4 = 0;b_step5 = 0;b_step6 = 0;b_tri = 0;b_filtered = 0;
  b_met = 0;
  b_weight = 1; b_puweight = 0; b_puweight_up = 0; b_puweight_dn = 0; b_genweight = 0; b_beffweight = 0; b_trigeffweight = 0;
  b_mueffweight = 0;b_mueffweight_up = 0;b_mueffweight_dn = 0;
  b_eleffweight = 0;b_eleffweight_up = 0;b_eleffweight_dn = 0;
  b_pdfWeights.clear();
  b_scaleWeights.clear();
  
  b_cat = 0; b_cat_eta = 0;
  b_isLoose = 0; b_isMedium = 0; b_isTight = 0;

  b_bumpcat1 = 0; b_bumpcat2 = 0;
 
  b_genlep1 = TLorentzVector(); b_genlep1_pid = 0;
  b_genlep2 = TLorentzVector(); b_genlep2_pid = 0;
  b_gendilep = TLorentzVector();
 
  b_lep1 = TLorentzVector(); b_lep1_pid = 0;
  b_lep2 = TLorentzVector(); b_lep2_pid = 0;
  b_jet1 = TLorentzVector(); b_jet2 = TLorentzVector();
  b_dilep = TLorentzVector(); b_dijet = TLorentzVector();
}

bool h2muAnalyzer::preSelection(const JetCollection& seljets, float met) const
{
  int njet = seljets.size();
  if (njet < 2) return false;

  if (seljets[0].pt() < 40.0) return false;
  
  if (seljets[1].pt() < 30.0) return false;
    
  //if (met > 40) return false;
  return true;
}

int h2muAnalyzer::jetCategory(const JetCollection& seljets, float met, float ll_pt) const
{  
  if (preSelection(seljets, met)){// pass preselection
    ofstream eventl;
    eventl.open("UOS_Event_List.txt", ios::app);
    eventl << b_event << endl;
    eventl.close();
    int a = seljets.size();
    //VBF Tight selection
    for (auto& j1 : seljets) {
        ofstream z;
        z.open("preselpass_dump.txt", ios::app);
        for (int i = 0; i < a; i++)
            {
            z<< "check jets1 "<<b_event<< " jet"<<i<<"->pt():" << seljets[i].pt() << " jet"<<i<<"->eta():" << seljets[i].eta() <<" jet"<<i<<"->phi():" << seljets[i].phi()<<endl;
            }
        z<<""<<endl;
        z.close();
        for (auto& j2 : seljets) {
	        if (&j1 != &j2) {//&& (j1.tlv().Pt() > 40.0) && (j2.tlv().Pt() > 30.0)){
	            TLorentzVector dijets = j1.tlv() + j2.tlv();
	            double delta_eta =  abs(j1.eta() - j2.eta());
                
	            if ((dijets.M() > 650 && abs(delta_eta) > 3.5)&& (j1.pt() >= 40 || j2.pt() >= 40)) {
                    ofstream file; 
                    file.open("UOS_VBF_Tight.txt", ios::app);
                    file << b_event<< endl;
                    file.close();
                    //int j = seljets.size();
                    ofstream VBFT;
                    VBFT.open("vbf-Event_Dump.txt", ios::app);
                    //VBFT<<"------------------------------------------------------------------"<<endl;
                    //VBFT<<b_event<<endl;
                    //VBFT<<""<<endl;
                    //for (int f = 0; f < j; f++){
                    //    VBFT<< "check jets1 "<<b_event<< " jet"<<f<<"->pt():" << seljets[f].pt() << " jet"<<f<<"->eta():" << seljets[f].eta() <<" jet"<<f<<"->phi():" << seljets[f].phi()<<endl;
                    //}
                    //VBFT<<""<<endl;
                    //VBFT<<"VBF_Tight "<<b_event<<" j1 pT:"<<j1.pt()<<" j2 pT:"<<j2.pt()<<" massJJ:"<<dijets.M()<<" DeltaJJ:"<<delta_eta<<endl;
                    VBFT<<b_event<<endl;
                    //VBFT<< " "<<endl;
                    VBFT.close();
  	                return 1; // VBF Tight selection
                }
            }
        }
    }
    for (auto& j3 : seljets) {
        for (auto& j4 : seljets) {
            if (&j3 !=  &j4){ 
                TLorentzVector dijets = j3.tlv() + j4.tlv();
                if ((dijets.M() > 250. && ll_pt >= 50.)&& (j3.pt() >= 40 || j4.pt() >= 40)){
                    ofstream file; 
                    file.open("UOS_ggF_Tight.txt", ios::app);
                    file << b_event<< endl;
                    file.close(); 
                    return 2; // GGF Tight selection 
                }
        	}
        }
    }
    ofstream file; 
    file.open("UOS_VBF_Loose.txt", ios::app);
    file << b_event<< endl;
    file.close(); 
    return 3; // VBF Loose selection    
  }
  else {
    ofstream eventl;
    eventl.open("UOS_Event_List.txt", ios::app);
    eventl << b_event << endl;
    eventl.close();
    if (ll_pt >= 25.){
      ofstream file; 
      file.open("UOS_01Jets_Tight.txt", ios::app);
      file << b_event<< endl;
      file.close(); 
      return 4; // 0,1-jet Tight
    }
    else {
      ofstream file; 
      file.open("UOS_01Jets_Loose.txt", ios::app);
      file << b_event<< endl;
      file.close(); 
      return 5; // 0,1-jet Loose
    }
  }
  
  return -1;
}

int h2muAnalyzer::etaCategory(float lep1_eta, float lep2_eta) const
{
  const float eta_mu[2] = {abs(lep1_eta),abs(lep2_eta)};
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
  
  return -1;
}

bool h2muAnalyzer::bumpCategory(const JetCollection& jets,const MuonCollection& recolep)
{
  /*
    for(int i=0;i<int(jets.size());i++){
    vector<double> row;
    row.push_back(jets[i].pt());
    row.push_back(abs(jets[i].eta()));
    row.push_back(jets[i].phi());
    row.push_back(jets[i].bDiscriminator(BTAG_CSVv2));
    jet2d.push_back(row);
    }
    sort(jet2d.begin(),jet2d.end(), [](vector<double> a, vector<double> b){return a[0] > b[0];});//sort by jet pt()
  */
  if (recolep[0].charge() * recolep[1].charge() > 0) return false;
  if (recolep[0].pt() < 25 || recolep[1].pt() < 25) return false;
  if (recolep[0].eta() > 2.1 || recolep[1].eta() > 2.1) return false;
  return true;
} 

bool h2muAnalyzer::bumpCat1(const JetCollection& jets,const MuonCollection& recolep)
{
  bool event_cat=false;
  if (!bumpCategory(jets,recolep)){return event_cat;}
  int btagpass=0;
  for(auto& j:jets){
    Jet jet(j);
    if (jet.bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2T) {
      //cout<<(*a)[0]<<", "<<(*a)[1]<<", "<<(*a)[2]<<", "<<WP_BTAG_CSVv2M<<endl;
      if ( jet.pt() > 30. && jet.eta() < 2.4 ) btagpass++;
    }
    else {
      if( jet.pt() > 30. && jet.eta() < 2.4 ) return false;
      if( jet.pt() > 30. && jet.eta() > 2.4 && jet.eta() < 4.7){
	event_cat=true;
      }
    }
  }
  if (btagpass==0){
    return false;
  }
  else {
    return event_cat;
  }
}

bool h2muAnalyzer::bumpCat2(const JetCollection& jets,const MuonCollection& recolep, float met)
{
  bool event_cat=false;
  if (!bumpCategory(jets,recolep)){return event_cat;}
  int btag=0;
  int cond=0;
  TLorentzVector dilep = recolep[0].tlv() + recolep[1].tlv();
  TLorentzVector dijet(0,0,0,0);
  for(auto& j:jets){
    Jet jet(j);
    if ( jet.pt() > 30. && jet.eta() > 2.4 && jet.eta() < 4.7 ) continue;
    if ( jet.pt() < 30. || jet.eta() > 2.4 ) continue;
    cond+=1;
    if ( jet.bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2T ){
      btag+=1;
      dijet+=jet.tlv();
    }
    else if ( jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2T ){
      dijet+=jet.tlv();
    }
  }
  //cout<<dijet.Phi()<<endl;
  // dilep and dijet passed the previous cut
  double delta_phi = abs(dilep.Phi()-dijet.Phi());
  // 2nd event category
  if ( (btag > 0) && (cond == 2) ){
    //printf("no error\n");
    if (met<40 && delta_phi>2.5){
      event_cat=true;
      //cout<<delta_phi<<endl;
    }
  }
  return event_cat; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(h2muAnalyzer);
