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
#include "CATTools/CatAnalyzer/src/rochcor2016.h"
#include "CATTools/CatAnalyzer/src/RoccoR.h"
//#include "TLorentzVector.h"

using namespace std;
using namespace cat;

class h2muAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit h2muAnalyzer(const edm::ParameterSet&);
  ~h2muAnalyzer();

  enum sys_e {sys_nom,
	      sys_jes_u, sys_jes_d, sys_jer_u, sys_jer_d,
	      sys_mu_u, sys_mu_d, sys_el_u, sys_el_d,
	      //sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
	      //sys_btag_u, sys_btag_d,
	      nsys_e
  };

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  void resetBr();
  float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, sys_e sys) const;
  float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, cat::MuonCollection& recolep, sys_e sys);
  int preSelect(const cat::JetCollection& seljets, float MET) const;
  int jetCategory(const cat::JetCollection& seljets, float MET, float ll_pt) const;
  int etaCategory(float lep1_eta, float lep2_eta) const;

  ScaleFactorEvaluator muonSF_, elecSF_;
  float getMuEffSF(const cat::Lepton& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 13 ) {
      const double pt = p.pt(), aeta = std::abs(p.eta());
      if      ( sys == +1 ) return muonSF_(pt, aeta,  1);
      else if ( sys == -1 ) return muonSF_(pt, aeta, -1);
      else return muonSF_(pt, aeta, 0);
    }
    return 1;
  }
  float getElEffSF(const cat::Lepton& p, int sys) const
  { 
    const int aid = abs(p.pdgId());
    if ( aid == 11 ) {
      const auto& el = dynamic_cast<const cat::Electron&>(p);
      const double pt = p.pt(), aeta = std::abs(el.scEta());
      if      ( sys == +1 ) return elecSF_(pt, aeta,  1);
      else if ( sys == -1 ) return elecSF_(pt, aeta, -1);
      else return elecSF_(pt, aeta, 0);
    }
    return 1;
  }
  
  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genweightToken_, puweightToken_, puweightToken_up_, puweightToken_dn_, topPtWeight_;
  //edm::EDGetTokenT<vector<float>> pdfweightToken_, scaleweightToken_;
  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;  
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  rochcor2016 *rmcor;
  bool runOnMC;
  
  std::vector<TTree*> ttree_;
  TH1D * h_nevents;
  int b_run, b_lumi, b_event;
  int b_nvertex, b_step, b_channel, b_njet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_filtered;
  float b_tri;
  float b_met, b_weight, b_puweight, b_puweight_up, b_puweight_dn, b_genweight,
    b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
    b_eleffweight, b_eleffweight_up, b_eleffweight_dn;
  std::vector<float> b_pdfWeights, b_scaleWeights;

  int b_cat;
  int b_cat_eta;
  bool b_isLoose, b_isMedium, b_isTight;

  TLorentzVector b_genlep1; int b_genlep1_pid;
  TLorentzVector b_genlep2; int b_genlep2_pid;
  TLorentzVector b_gendilep;
  
  TLorentzVector b_lep1; int b_lep1_pid;
  TLorentzVector b_lep2; int b_lep2_pid;
  TLorentzVector b_jet1, b_jet2, b_dilep, b_dijet;
  
  const static int NCutflow = 7;
  std::vector<std::vector<int> > cutflow_;
};
//
// constructors and destructor
//
h2muAnalyzer::h2muAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  //genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  //pdfweightToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("pdfweight"));
  //scaleweightToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaleweight"));
  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightToken_up_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_up"));
  puweightToken_dn_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_dn"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  typedef std::vector<double> vdouble;

  const auto muonSet = iConfig.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("effSF");
  muonSF_.set(muonSFSet.getParameter<vdouble>("pt_bins"),
              muonSFSet.getParameter<vdouble>("abseta_bins"),
              muonSFSet.getParameter<vdouble>("values"),
              muonSFSet.getParameter<vdouble>("errors"));

  const auto elecSet = iConfig.getParameter<edm::ParameterSet>("electron");
  elecToken_ = consumes<cat::ElectronCollection>(elecSet.getParameter<edm::InputTag>("src"));
  const auto elecSFSet = elecSet.getParameter<edm::ParameterSet>("effSF");
  elecSF_.set(elecSFSet.getParameter<vdouble>("pt_bins"),
              elecSFSet.getParameter<vdouble>("abseta_bins"),
              elecSFSet.getParameter<vdouble>("values"),
              elecSFSet.getParameter<vdouble>("errors"));


  usesResource("TFileService");
  edm::Service<TFileService> fs;
  const std::string sys_name[nsys_e] = {
    "nom",
    "jes_u", "jes_d", "jer_u", "jer_d",
    "mu_u", "mu_d", "el_u", "el_d",
    //    "mueff_u", "mueff_d", "eleff_u", "eleff_d",
    //    "btag_u", "btag_d"
  };

  h_nevents = fs->make<TH1D>("nevents","nevents",1,0,1);       
  for (int sys = 0; sys < nsys_e; ++sys){
    ttree_.push_back(fs->make<TTree>(sys_name[sys].c_str(), sys_name[sys].c_str()));
    auto tr = ttree_.back();
    tr->Branch("run", &b_run, "run/I");
    tr->Branch("event", &b_event, "event/I");

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
    
    if (sys == 0){
      tr->Branch("weight", &b_weight, "weight/F");
      tr->Branch("puweight", &b_puweight, "puweight/F");
      tr->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
      tr->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
      tr->Branch("genweight", &b_genweight, "genweight/F");
      tr->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
      tr->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
      tr->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
      tr->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
      tr->Branch("eleffweight_up", &b_eleffweight_up, "eleffweight_up/F");
      tr->Branch("eleffweight_dn", &b_eleffweight_dn, "eleffweight_dn/F");      
      tr->Branch("pdfWeights","std::vector<float>",&b_pdfWeights);
      tr->Branch("scaleWeights","std::vector<float>",&b_scaleWeights);
      
      tr->Branch("genlep1", "TLorentzVector", &b_genlep1);
      tr->Branch("genlep1_pid", &b_genlep1_pid, "genlep1_pid/I");    
      tr->Branch("genlep2", "TLorentzVector", &b_genlep2);
      tr->Branch("genlep2_pid", &b_genlep2_pid, "genlep2_pid/I");    
      tr->Branch("gendilep", "TLorentzVector", &b_gendilep);
    }
  }
 
  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});
}

h2muAnalyzer::~h2muAnalyzer()
{
  cout <<"     cut flow   emu    ee    mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step"<< i << "    "<< cutflow_[i][0] <<  "   "<< cutflow_[i][1] << "   " << cutflow_[i][2] << "   " << cutflow_[i][3]<< endl;
  }
}

void h2muAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&){}

void h2muAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  rmcor = new rochcor2016();  
  runOnMC = !iEvent.isRealData();
  
  cutflow_[0][0]++;

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();
    
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

      // edm::Handle<float> genweightHandle;
      // iEvent.getByToken(genweightToken_, genweightHandle);
      // b_genweight = (*genweightHandle);
      b_weight = b_genweight*b_puweight;

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
    
    if (sys == sys_nom)
      h_nevents->Fill(0.5,b_puweight*b_genweight);

    edm::Handle<int> lumiSelectionHandle;
    iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
    if (!runOnMC){
      if (*lumiSelectionHandle == 0) return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    // filters
    edm::Handle<int> recoFiltersHandle;
    iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
    b_filtered = *recoFiltersHandle == 0 ? false : true;
    if (b_filtered){
      b_step1 = true;
      b_step = 1;
      if (sys == sys_nom) cutflow_[1][b_channel]++;
    }
    ////////////////////////////////////////////////////////////////////////////////
    // trigger
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
    AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

    if (trigHelper.triggerFired("HLT_IsoMu20_v") || trigHelper.triggerFired("HLT_IsoTrkMu20_v")){
      b_step2 = true;
      if (b_step == 1){
	b_step = 2;
	if (sys == sys_nom) cutflow_[2][b_channel]++;
      }
    }
    ////////////////////////////////////////////////////////////////////////////////
    // vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (!vertices->empty()){ // skip the event if no PV found
      b_step3 = true;
      if (b_step == 2){
	b_step = 3;
	if (sys == sys_nom) cutflow_[3][b_channel]++;
      }
    }
    // const reco::Vertex &PV = vertices->front();
    edm::Handle<int> nGoodVertexHandle;
    iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
    b_nvertex = *nGoodVertexHandle;
    
    ////////////////////////////////////////////////////////////////////////////////
    // muon selection
    edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
    cat::MuonCollection selectedMuons;
    selectMuons(*muons, selectedMuons, (sys_e)sys);
    if ( selectedMuons.size() < 2 ) continue;
    if (selectedMuons[0].charge()*selectedMuons[1].charge() > 0) continue;

    b_lep1 = selectedMuons[0].tlv();b_lep1_pid = selectedMuons[0].pdgId();
    b_lep2 = selectedMuons[1].tlv();b_lep2_pid = selectedMuons[1].pdgId();
    b_dilep = b_lep1 + b_lep2;
    
    b_isLoose = (selectedMuons[0].isLooseMuon() && selectedMuons[1].isLooseMuon());
    b_isMedium = (selectedMuons[0].isMediumMuon() && selectedMuons[1].isMediumMuon());
    b_isTight = (selectedMuons[0].isTightMuon() && selectedMuons[1].isTightMuon());
    
    if (b_isTight){
      b_step4 = true;
      if (b_step == 3){
	b_step = 4;
	if (sys == sys_nom) cutflow_[4][b_channel]++;
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // trigger matching
    bool triggerMatching = false;
    if ( selectedMuons[0].pt() > 20 &&
	 (trigHelper.triggerMatched("HLT_IsoMu20_v", selectedMuons[0]) ||
	  trigHelper.triggerMatched("HLT_IsoTrkMu20_v", selectedMuons[0])))
      triggerMatching = true;
    
    if ( selectedMuons[1].pt() > 20 &&
	 (trigHelper.triggerMatched("HLT_IsoMu20_v", selectedMuons[1]) ||
	  trigHelper.triggerMatched("HLT_IsoTrkMu20_v", selectedMuons[1])))
      triggerMatching = true;

    if (triggerMatching){
      b_tri = true;
      b_step5 = true;
      if (b_step == 4){
	b_step = 5;
	if (sys == sys_nom) cutflow_[5][b_channel]++;
      }
    }
    ////////////////////////////////////////////////////////////////////////////////
    // jets
    edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
    JetCollection&& selectedJets = selectJets(*jets, selectedMuons, (sys_e)sys);

    if (selectedJets.size() <= 2){
      b_step6 = true;
      if (b_step == 5){
	b_step = 6;
	if (sys == sys_nom) cutflow_[6][b_channel]++;
      }
    }
    
    b_mueffweight    = getMuEffSF(selectedMuons[0],  0)*getMuEffSF(selectedMuons[1],  0);
    b_mueffweight_up = getMuEffSF(selectedMuons[0], +1)*getMuEffSF(selectedMuons[1], +1);
    b_mueffweight_dn = getMuEffSF(selectedMuons[0], -1)*getMuEffSF(selectedMuons[1], -1);
    
    b_njet = selectedJets.size();
    if (b_njet > 0)
      b_jet1 = selectedJets[0].tlv();
    if (b_njet > 1){
      b_jet2 = selectedJets[1].tlv();
      b_dijet = b_jet1 + b_jet2;
    }

    edm::Handle<cat::METCollection> mets; iEvent.getByToken(metToken_, mets);
    const auto met = mets->front().p4();
    b_met = met.pt();

    b_cat_eta = etaCategory(b_lep1.Eta(), b_lep2.Eta());
    b_cat = jetCategory(selectedJets, b_met, b_dilep.Pt());

    ttree_[sys]->Fill();
  }
  delete rmcor;
}

float h2muAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, sys_e sys) const
{
  float weight = 1.;
  for (auto& m : muons) {
    cat::Muon mu(m);
    if (std::abs(mu.eta()) > 2.4) continue;
    
    float qter = 1.0;
    TLorentzVector tmu(mu.tlv());

    cout << "\n initial   " << mu.tlv().Pt()
	 << ", " << mu.tlv().Eta()
	 << ", " << mu.tlv().Phi()
	 << ", " << mu.tlv().M()
	 << " mu.trackerLayersWithMeasurement() " << mu.trackerLayersWithMeasurement()
	 << " mu.numberOfValidHits() " << mu.numberOfValidHits()
	 <<endl;
    
    if (runOnMC)
      rmcor->momcor_mc(tmu, mu.charge(), mu.trackerLayersWithMeasurement(), qter);
    else 
      rmcor->momcor_data(tmu, mu.charge(), 0, qter);

    mu.setP4(m.p4() * tmu.E()/mu.tlv().E());

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
    if (sys == sys_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == sys_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());

    if (mu.pt() < 10.) continue;
    if (!mu.isLooseMuon()) continue;
    if (mu.relIso(0.4) > 0.25) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    //weight *= mu.scaleFactor("NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1");
    //weight *= mu.scaleFactor("NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1");
    selmuons.push_back(mu);
  }
  return weight;
}

float h2muAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
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

cat::JetCollection h2muAnalyzer::selectJets(const cat::JetCollection& jets, cat::MuonCollection& recolep, sys_e sys)
{
  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    if (sys == sys_jes_u) jet.setP4(j.p4() * j.shiftedEnUp());
    if (sys == sys_jes_d) jet.setP4(j.p4() * j.shiftedEnDown());
    if (sys == sys_jer_u) jet.setP4(j.p4() * j.smearedResUp());
    if (sys == sys_jer_d) jet.setP4(j.p4() * j.smearedResDown());

    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 4.7)  continue;
    //if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep.p4()) < 0.3) hasOverLap = true;
    }
    if (hasOverLap) continue;
    seljets.push_back(jet);
  }

  return seljets;
}

void h2muAnalyzer::resetBr()
{
  b_nvertex = 0;b_step = 0;b_channel = 0;b_njet = 0;
  b_step1 = 0;b_step2 = 0;b_step3 = 0;b_step4 = 0;b_step5 = 0;b_step6 = 0;b_tri = 0;b_filtered = 0;
  b_met = 0;
  b_weight = 0; b_puweight = 0; b_puweight_up = 0; b_puweight_dn = 0; b_genweight = 0;
  b_mueffweight = 0;b_mueffweight_up = 0;b_mueffweight_dn = 0;
  b_eleffweight = 0;b_eleffweight_up = 0;b_eleffweight_dn = 0;
  b_pdfWeights.clear();
  b_scaleWeights.clear();

  b_cat = 0; b_cat_eta = 0;
  b_isLoose = 0; b_isMedium = 0; b_isTight = 0;
 
  b_genlep1 = TLorentzVector(); b_genlep1_pid = 0;
  b_genlep2 = TLorentzVector(); b_genlep2_pid = 0;
  b_gendilep = TLorentzVector();

  b_lep1 = TLorentzVector(); b_lep1_pid = 0;
  b_lep2 = TLorentzVector(); b_lep2_pid = 0;
  b_jet1 = TLorentzVector(); b_jet2 = TLorentzVector();
  b_dilep = TLorentzVector(); b_dijet = TLorentzVector();
}
int h2muAnalyzer::preSelect(const cat::JetCollection& seljets, float met) const
{
  int njet = seljets.size();
  if (njet == 2){
    if (seljets[0].pt()>40 && seljets[1].pt()>30 && met<40){
      return 1; // pass preselection
    }
  }
  return 0;// fail preselection
}

int h2muAnalyzer::jetCategory(const cat::JetCollection& seljets, float met, float ll_pt) const
{
  int presel = preSelect(seljets, met);
  
  if (presel == 0){// fail preselection
    if (ll_pt >= 10) return 4; // 0,1-jet Tight
    else             return 5; // 0,1-jet Loose
  }
  
  if (presel == 1){// pass preselection
    TLorentzVector M_jets = seljets[0].tlv() + seljets[1].tlv();
    double delta_eta = seljets[0].eta()-seljets[1].eta();

    bool VBF_Tight = (M_jets.M() > 650 && abs(delta_eta) > 3.5);
    bool ggF_Tight = (M_jets.M() > 250 && ll_pt > 50);

    if (VBF_Tight) return 1; // VBF Tight selection
    else {
      if (ggF_Tight) return 2; // ggF Tight selection
      else return 3; // VBF Loose selection
    }
  } 
  return -1;
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
  
  return -1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(h2muAnalyzer);
