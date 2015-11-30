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

//#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"
#include "CATTools/CatAnalyzer/interface/CSVHelper.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"
//#include "TLorentzVector.h"

using namespace std;
using namespace cat;

class TtbarBbbarDiLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit TtbarBbbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarBbbarDiLeptonAnalyzer();

  enum {
    CH_NONE=0, CH_MUEL=1, CH_ELEL=2, CH_MUMU=3
  };

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  void selectMuons(const cat::MuonCollection& muons, LeptonCollection& selmuons) const;
  void selectElecs(const cat::ElectronCollection& elecs, LeptonCollection& selelecs) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonCollection& recolep) const;
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  const reco::Candidate* getLast(const reco::Candidate* p) const;

  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_;
  edm::EDGetTokenT<float> genweightToken_, puweightToken_;
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<vector<int> > partonTop_modes_;

  edm::EDGetTokenT<vector<int> > genTtbarIdToken_, NgenJet30Token_, genTtbarLeptonDecayToken_;

  edm::EDGetTokenT<reco::GenParticleCollection> partonTop_genParticles_, pseudoTop_;

  TTree * ttree_;
  int b_nvertex, b_step, b_channel;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_tri, b_filtered;
  float b_met, b_metphi;
  int b_njet30, b_nbjetM30, b_nbjetT30;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi, b_lep1_RelIso;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi, b_lep2_RelIso;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;
 
 //for ttbb 
  std::vector<double> b_jets_pt;
  std::vector<double> b_jets_eta;
  std::vector<double> b_jets_phi;
  std::vector<int>    b_jets_flavor;
  std::vector<double> b_jets_bDiscriminatorCSV;
  std::vector<int>    b_csvd_jetid;
 
  //mc
  float b_weight, b_puweight;
  int b_partonChannel, b_partonMode1, b_partonMode2;
  float b_partonlep1_pt, b_partonlep1_eta;
  float b_partonlep2_pt, b_partonlep2_eta;
  bool b_partonInPhase, b_partonInPhaseJet, b_partonInPhaseLep;
  int b_pseudoTopChannel;
  float b_pseudoToplep1_pt, b_pseudoToplep1_eta;
  float b_pseudoToplep2_pt, b_pseudoToplep2_eta;
  bool b_pseudoInPhase;

  int b_genTtbarId, b_NgenJet30, b_genTtbarLeptonDecay;

  //mc
  float  b_csvweight;
  float  b_csvweight_JES_Up;
  float  b_csvweight_JES_Down;
  float  b_csvweight_LF_Up;
  float  b_csvweight_LF_Down;
  float  b_csvweight_HF_Up;
  float  b_csvweight_HF_Down;
  float  b_csvweight_HF_Stats1_Up;
  float  b_csvweight_HF_Stats1_Down;
  float  b_csvweight_HF_Stats2_Up;
  float  b_csvweight_HF_Stats2_Down;
  float  b_csvweight_LF_Stats1_Up;
  float  b_csvweight_LF_Stats1_Down;
  float  b_csvweight_LF_Stats2_Up;
  float  b_csvweight_LF_Stats2_Down;
  float  b_csvweight_Charm_Err1_Up;
  float  b_csvweight_Charm_Err1_Down;
  float  b_csvweight_Charm_Err2_Up;
  float  b_csvweight_Charm_Err2_Down;


/* 
  //float b_jet1_pt, b_jet1_eta, b_jet1_CSVInclV2;
  //float b_jet2_pt, b_jet2_eta, b_jet2_CSVInclV2;
  float b_top1_pt, b_top1_eta, b_top1_phi, b_top1_rapi;
  float b_top2_pt, b_top2_eta, b_top2_phi, b_top2_rapi;
  float b_ttbar_pt, b_ttbar_eta, b_ttbar_phi, b_ttbar_m, b_ttbar_rapi;
  float b_maxweight;*/
  int b_is3lep;

  //std::unique_ptr<TtFullLepKinSolver> solver;
  std::unique_ptr<KinematicSolver> solver_;
  //enum TTbarMode { CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON };
  //enum DecayMode { CH_HADRON = 1, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

  const static int NCutflow = 10;
  std::vector<std::vector<int> > cutflow_;
  bool runOnMC_;
  CSVHelper *csvWeight_;
};
//
// constructors and destructor
//
TtbarBbbarDiLeptonAnalyzer::TtbarBbbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  trigTokenMUEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));

//  genTtbarIdToken_, NgenJet30Token_, genTtbarLeptonDecayToken_;
  genTtbarIdToken_          = consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarId"));
  NgenJet30Token_           = consumes<int>(iConfig.getParameter<edm::InputTag>("NgenJet30"));
  genTtbarLeptonDecayToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarLeptonDecay"));

  muonToken_ = consumes<cat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  partonTop_channel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("partonTop_channel"));
  partonTop_modes_   = consumes<vector<int> >(iConfig.getParameter<edm::InputTag>("partonTop_modes"));
  partonTop_genParticles_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("partonTop_genParticles"));
  pseudoTop_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pseudoTop"));

  csvWeight = new CSVHelper();
/*
  auto solverPSet = iConfig.getParameter<edm::ParameterSet>("solver");
  auto algoName = solverPSet.getParameter<std::string>("algo");
  std::transform(algoName.begin(), algoName.end(), algoName.begin(), ::toupper);
  if      ( algoName == "CMSKIN" ) solver_.reset(new CMSKinSolver(solverPSet));
  else if ( algoName == "DESYMASSLOOP" ) solver_.reset(new DESYMassLoopSolver(solverPSet));
  else if ( algoName == "DESYSMEARED" ) solver_.reset(new DESYSmearedSolver(solverPSet));
  else if ( algoName == "MT2"    ) solver_.reset(new MT2Solver(solverPSet));
  else if ( algoName == "MAOS"   ) solver_.reset(new MAOSSolver(solverPSet));
  else if ( algoName == "NUWGT"  ) solver_.reset(new NuWeightSolver(solverPSet));
  else if ( algoName == "DEFAULT" ) solver_.reset(new TTDileptonSolver(solverPSet));
  else {
    cerr << "The solver name \"" << solverPSet.getParameter<std::string>("algo") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solver_.reset(new TTDileptonSolver(solverPSet)); // A dummy solver
  }
*/

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("nom", "nom");
  ttree_->Branch("nvertex", &b_nvertex, "nvertex/I");
  ttree_->Branch("step", &b_step, "step/I");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("njet30", &b_njet30, "njet30/I");
  ttree_->Branch("nbjetM30", &b_nbjetM30, "nbjetM30/I");
  ttree_->Branch("nbjetT30", &b_nbjetT30, "nbjetT30/I");
  ttree_->Branch("step1", &b_step1, "step1/O");
  ttree_->Branch("step2", &b_step2, "step2/O");
  ttree_->Branch("step3", &b_step3, "step3/O");
  ttree_->Branch("step4", &b_step4, "step4/O");
  ttree_->Branch("step5", &b_step5, "step5/O");
  ttree_->Branch("step6", &b_step6, "step6/O");
  ttree_->Branch("tri", &b_tri, "tri/O");
  ttree_->Branch("filtered", &b_filtered, "filtered/O");
  ttree_->Branch("met", &b_met, "met/F");
  ttree_->Branch("metphi", &b_metphi, "metphi/F");

  ttree_->Branch("weight", &b_weight, "weight/F");
  ttree_->Branch("puweight", &b_puweight, "puweight/F");

  ttree_->Branch("lep1_pt", &b_lep1_pt, "lep1_pt/F");
  ttree_->Branch("lep1_eta", &b_lep1_eta, "lep1_eta/F");
  ttree_->Branch("lep1_phi", &b_lep1_phi, "lep1_phi/F");
  ttree_->Branch("lep1_RelIso", &b_lep1_RelIso, "lep1_RelIso/F");
  ttree_->Branch("lep2_pt", &b_lep2_pt, "lep2_pt/F");
  ttree_->Branch("lep2_eta", &b_lep2_eta, "lep2_eta/F");
  ttree_->Branch("lep2_phi", &b_lep2_phi, "lep2_phi/F");
  ttree_->Branch("lep2_RelIso", &b_lep2_RelIso, "lep2_RelIso/F");
  ttree_->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
  ttree_->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  ttree_->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  ttree_->Branch("ll_m", &b_ll_m, "ll_m/F");
  
  ttree_->Branch("parton_channel", &b_partonChannel, "parton_channel/I");
  ttree_->Branch("parton_mode1", &b_partonMode1, "parton_mode1/I");
  ttree_->Branch("partonlep1_pt", &b_partonlep1_pt, "partonlep1_pt/F");
  ttree_->Branch("partonlep1_eta", &b_partonlep1_eta, "partonlep1_eta/F");
  ttree_->Branch("partonlep2_pt", &b_partonlep2_pt, "partonlep2_pt/F");
  ttree_->Branch("partonlep2_eta", &b_partonlep2_eta, "partonlep2_eta/F");
  ttree_->Branch("parton_mode2", &b_partonMode2, "parton_mode2/I");
  ttree_->Branch("partonInPhase", &b_partonInPhase, "partonInPhase/O");
  ttree_->Branch("partonInPhaseLep", &b_partonInPhaseLep, "partonInPhaseLep/O");
  ttree_->Branch("partonInPhaseJet", &b_partonInPhaseJet, "partonInPhaseJet/O");

  ttree_->Branch("pseudoTop_channel", &b_pseudoTopChannel, "pseudoTop_channel/I");
  ttree_->Branch("pseudoToplep1_pt", &b_pseudoToplep1_pt, "pseudoToplep1_pt/F");
  ttree_->Branch("pseudoToplep1_eta", &b_pseudoToplep1_eta, "pseudoToplep1_eta/F");
  ttree_->Branch("pseudoToplep2_pt", &b_pseudoToplep2_pt, "pseudoToplep2_pt/F");
  ttree_->Branch("pseudoToplep2_eta", &b_pseudoToplep2_eta, "pseudoToplep2_eta/F");
  ttree_->Branch("pseudoInPhase", &b_pseudoInPhase, "pseudoInPhase/O");

  ttree->Branch("jets_pt","std::vector<double>",&b_jets_pt);
  ttree->Branch("jets_eta","std::vector<double>",&b_jets_eta);
  ttree->Branch("jets_phi","std::vector<double>",&b_jets_phi);
  ttree->Branch("jets_flavor","std::vector<int>",&b_jets_flavor);
  ttree->Branch("jets_bDiscriminatorCSV","std::vector<double>",&b_jets_bDiscriminatorCSV);
  ttree->Branch("csvd_jetid","std::vector<int>",&b_csvd_jetid);

  ttree_->Branch("csvweight", &b_csvweight, "csvweight/F");
  ttree_->Branch("csvweight_JES_Up",          &b_csvweight_JES_Up,          "csvweight_JES_Up/F");          
  ttree_->Branch("csvweight_JES_Down",        &b_csvweight_JES_Down,        "csvweight_JES_Down/F");
  ttree_->Branch("csvweight_LF_Up",           &b_csvweight_LF_Up,           "csvweight_LF_Up/F");
  ttree_->Branch("csvweight_LF_Down",         &b_csvweight_LF_Down,         "csvweight_LF_Down/F");
  ttree_->Branch("csvweight_HF_Up",           &b_csvweight_HF_Up,           "csvweight_HF_Up/F");
  ttree_->Branch("csvweight_HF_Down",         &b_csvweight_HF_Down,         "csvweight_HF_Down/F");
  ttree_->Branch("csvweight_HF_Stats1_Up",    &b_csvweight_HF_Stats1_Up,    "csvweight_HF_Stats1_Up/F");
  ttree_->Branch("csvweight_HF_Stats1_Down",  &b_csvweight_HF_Stats1_Down,  "csvweight_HF_Stats1_Down/F");
  ttree_->Branch("csvweight_HF_Stats2_Up",    &b_csvweight_HF_Stats2_Up,    "csvweight_HF_Stats2_Up/F");
  ttree_->Branch("csvweight_HF_Stats2_Down",  &b_csvweight_HF_Stats2_Down,  "csvweight_HF_Stats2_Down/F");
  ttree_->Branch("csvweight_LF_Stats1_Up",    &b_csvweight_LF_Stats1_Up,    "csvweight_LF_Stats1_Up/F");
  ttree_->Branch("csvweight_LF_Stats1_Down",  &b_csvweight_LF_Stats1_Down,  "csvweight_LF_Stats1_Down/F");
  ttree_->Branch("csvweight_LF_Stats2_Up",    &b_csvweight_LF_Stats2_Up,    "csvweight_LF_Stats2_Up/F");
  ttree_->Branch("csvweight_LF_Stats2_Down",  &b_csvweight_LF_Stats2_Down,  "csvweight_LF_Stats2_Down/F");
  ttree_->Branch("csvweight_Charm_Err1_Up",   &b_csvweight_Charm_Err1_Up,   "csvweight_Charm_Err1_Up/F");
  ttree_->Branch("csvweight_Charm_Err1_Down", &b_csvweight_Charm_Err1_Down, "csvweight_Charm_Err1_Down/F");
  ttree_->Branch("csvweight_Charm_Err2_Up",   &b_csvweight_Charm_Err2_Up,   "csvweight_Charm_Err2_Up/F");
  ttree_->Branch("csvweight_Charm_Err2_Down", &b_csvweight_Charm_Err2_Down, "csvweight_Charm_Err2_Down/F");


/*
  ttree_->Branch("jet1_pt", &b_jet1_pt, "jet1_pt/F");
  ttree_->Branch("jet2_pt", &b_jet2_pt, "jet2_pt/F");
  ttree_->Branch("jet1_eta", &b_jet1_eta, "jet1_eta/F");
  ttree_->Branch("jet2_eta", &b_jet2_eta, "jet2_eta/F");
  ttree_->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
  ttree_->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");

  ttree_->Branch("top1_pt", &b_top1_pt, "top1_pt/F");
  ttree_->Branch("top1_eta", &b_top1_eta, "top1_eta/F");
  ttree_->Branch("top1_phi", &b_top1_phi, "top1_phi/F");
  ttree_->Branch("top1_rapi", &b_top1_rapi, "top1_rapi/F");
  ttree_->Branch("top2_pt", &b_top2_pt, "top2_pt/F");
  ttree_->Branch("top2_eta", &b_top2_eta, "top2_eta/F");
  ttree_->Branch("top2_phi", &b_top2_phi, "top2_phi/F");
  ttree_->Branch("top2_rapi", &b_top2_rapi, "top2_rapi/F");
  ttree_->Branch("ttbar_pt", &b_ttbar_pt, "ttbar_pt/F");
  ttree_->Branch("ttbar_eta", &b_ttbar_eta, "ttbar_eta/F");
  ttree_->Branch("ttbar_phi", &b_ttbar_phi, "ttbar_phi/F");
  ttree_->Branch("ttbar_rapi", &b_ttbar_rapi, "ttbar_rapi/F");
  ttree_->Branch("ttbar_m", &b_ttbar_m, "ttbar_m/F");
*/
  ttree_->Branch("is3lep", &b_is3lep, "is3lep/I");

  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});
}

TtbarBbbarDiLeptonAnalyzer::~TtbarBbbarDiLeptonAnalyzer()
{
  cout <<"cut flow         emu         ee         mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step "<< i << " "<< cutflow_[i][0] <<  " "<< cutflow_[i][1] << " " << cutflow_[i][2] << " " << cutflow_[i][3]<< endl;
  }
}

void TtbarBbbarDiLeptonAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&)
{
  if ( dynamic_cast<DESYSmearedSolver*>(solver_.get()) != 0 ) {
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(lumi.index());
    dynamic_cast<DESYSmearedSolver*>(solver_.get())->setRandom(&engine);
  }
}

void TtbarBbbarDiLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  b_nvertex = 0;b_step = -1;b_channel = 0;
  b_njet30 = 0; b_nbjetM30 = 0; b_nbjetT30 = 0;
  b_step1 = 0;b_step2 = 0;b_step3 = 0;b_step4 = 0;b_step5 = 0;b_step6 = 0;b_tri = 0;b_filtered = 0;
  b_met = -9; b_metphi = -9;
  b_weight = 1; b_puweight = 1;
  
  b_lep1_pt = -9;b_lep1_eta = -9;b_lep1_phi = -9; b_lep1_RelIso = -9;
  b_lep2_pt = -9;b_lep2_eta = -9;b_lep2_phi = -9; b_lep2_RelIso = -9;
  b_ll_pt = -9;b_ll_eta = -9;b_ll_phi = -9;b_ll_m = -9;

  b_partonChannel = -1; b_partonMode1 = -1; b_partonMode2 = -1;
  b_partonlep1_pt = -9; b_partonlep1_eta = -9;
  b_partonlep2_pt = -9; b_partonlep2_eta = -9;
  b_partonInPhase = 0; b_partonInPhaseLep = false; b_partonInPhaseJet = false;
  b_pseudoTopChannel = -1;
  b_pseudoToplep1_pt = -9; b_pseudoToplep1_eta = -9;
  b_pseudoToplep2_pt = -9; b_pseudoToplep2_eta = -9;
  b_pseudoInPhase = false;

  b_jets_pt->clear();
  b_jets_eta->clear();                 
  b_jets_phi->clear();
  b_jets_flavor->clear();
  b_jets_bDiscriminatorCSV->clear();
  b_csvd_jetid->clear();

  b_csvweight = 1;
  b_csvweight_JES_Up = 1;
  b_csvweight_JES_Down = 1;
  b_csvweight_LF_Up = 1;
  b_csvweight_LF_Down = 1;
  b_csvweight_HF_Up = 1;
  b_csvweight_HF_Down = 1;
  b_csvweight_HF_Stats1_Up = 1;
  b_csvweight_HF_Stats1_Down = 1;
  b_csvweight_HF_Stats2_Up = 1;
  b_csvweight_HF_Stats2_Down = 1;
  b_csvweight_LF_Stats1_Up = 1;
  b_csvweight_LF_Stats1_Down = 1;
  b_csvweight_LF_Stats2_Up = 1;
  b_csvweight_LF_Stats2_Down = 1;
  b_csvweight_Charm_Err1_Up = 1;
  b_csvweight_Charm_Err1_Down = 1;
  b_csvweight_Charm_Err2_Up = 1;
  b_csvweight_Charm_Err2_Down = 1;


/*  
  b_jet1_pt = -9; b_jet1_eta = -9; b_jet1_CSVInclV2 = -9;
  b_jet2_pt = -9; b_jet2_eta = -9; b_jet2_CSVInclV2 = -9;
  b_top1_pt = -9; b_top1_eta = -9; b_top1_phi = -9; b_top1_rapi = -9;
  b_top2_pt = -9; b_top2_eta = -9; b_top2_phi = -9; b_top2_rapi = -9;
  b_ttbar_pt = -9; b_ttbar_eta = -9; b_ttbar_phi = -9; b_ttbar_m = -9; b_ttbar_rapi = -9;
*/
  b_is3lep = -9;

  cutflow_[0][b_channel]++;
  
  edm::Handle<int> partonTop_channel;
  if ( iEvent.getByToken(partonTop_channel_, partonTop_channel)){
    edm::Handle<vector<int> > partonTop_modes;
    edm::Handle<reco::GenParticleCollection> partonTop_genParticles;
    iEvent.getByToken(partonTop_modes_, partonTop_modes);
    iEvent.getByToken(partonTop_genParticles_, partonTop_genParticles);
    if ( (*partonTop_modes).size() == 0 ) {
      b_partonMode1 = 0;
      b_partonMode2 = 0;
    }
    else if ( (*partonTop_modes).size() == 1 ) { b_partonMode2 = 0; }
    else{
      b_partonChannel = *partonTop_channel;
      b_partonMode1 = (*partonTop_modes)[0];
      b_partonMode2 = (*partonTop_modes)[1];
    }

    if ( !(partonTop_genParticles->empty()) ){

      // Get Top quark pairs
      const auto parton1 = &partonTop_genParticles->at(0);
      const auto parton2 = &partonTop_genParticles->at(1);
      // Get W and b quarks
      if ( parton1 and parton2 ) {
        const auto partonW1 = parton1->daughter(0);
        const auto partonB1 = parton1->daughter(1);
        const auto partonW2 = parton2->daughter(0);
        const auto partonB2 = parton2->daughter(1);

	if ( (partonB1->pt() > 30 && std::abs(partonB1->eta()) < 2.4) && 
	     (partonB2->pt() > 30 && std::abs(partonB2->eta()) < 2.4))
	  b_partonInPhaseJet = true;
	
        // Get W daughters
        if ( partonW1 and partonW2 and partonB1 and partonB2 ) {
          const auto partonW11 = partonW1->daughter(0);
          const auto partonW21 = partonW2->daughter(0);
	  if ( (partonW11->pt() > 20 && std::abs(partonW11->eta()) < 2.4 && (std::abs(partonW11->pdgId()) == 11 || std::abs(partonW11->pdgId()) == 13) ) && 
	       (partonW21->pt() > 20 && std::abs(partonW21->eta()) < 2.4 && (std::abs(partonW11->pdgId()) == 11 || std::abs(partonW11->pdgId()) == 13) ))
	    b_partonInPhaseLep = true;

          // Fill lepton informations
          b_partonlep1_pt = partonW11->pt();
          b_partonlep1_eta = partonW11->eta();
          b_partonlep2_pt = partonW21->pt();
          b_partonlep2_eta = partonW21->eta();
        }
      }
      if (b_partonInPhaseJet && b_partonInPhaseLep) b_partonInPhase = true;
    }

    edm::Handle<reco::GenParticleCollection> pseudoTopHandle;
    iEvent.getByToken(pseudoTop_          , pseudoTopHandle);
    if ( !(pseudoTopHandle->empty()) ){
      b_pseudoTopChannel = CH_NONE;

      // Get Top quark pairs
      const auto pseudoTop1 = &pseudoTopHandle->at(0);
      const auto pseudoTop2 = &pseudoTopHandle->at(1);

      // Get W and b quarks
      if ( pseudoTop1 and pseudoTop2 ) {
        const auto pseudoW1 = pseudoTop1->daughter(0);
        const auto pseudoB1 = pseudoTop1->daughter(1);
        const auto pseudoW2 = pseudoTop2->daughter(0);
        const auto pseudoB2 = pseudoTop2->daughter(1);

        // Get W daughters
        if ( pseudoW1 and pseudoW2 and pseudoB1 and pseudoB2 ) {
          const auto pseudoW11 = pseudoW1->daughter(0);
          const auto pseudoW21 = pseudoW2->daughter(0);

          // Fill leps informations
          const int pseudoW1DauId = abs(pseudoW11->pdgId());
          const int pseudoW2DauId = abs(pseudoW21->pdgId());
          b_pseudoToplep1_pt = pseudoW11->pt();
          b_pseudoToplep1_eta = pseudoW11->eta();
          b_pseudoToplep2_pt = pseudoW21->pt();
          b_pseudoToplep2_eta = pseudoW21->eta();
          if ( pseudoW1DauId > 10 and pseudoW2DauId > 10 ) {
            switch ( pseudoW1DauId+pseudoW2DauId ) {
	    case 22: b_pseudoTopChannel = CH_ELEL; break;
	    case 26: b_pseudoTopChannel = CH_MUMU; break;
	    default: b_pseudoTopChannel = CH_MUEL;
            }
          }
	  b_partonInPhase = true;
        }
      }
    }
  }

  if (runOnMC_){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    b_puweight = *puweightHandle;
    edm::Handle<float> genweightHandle;
    iEvent.getByToken(genweightToken_, genweightHandle);
    b_weight = (*genweightHandle)*b_puweight;

    //////
    edm::Handle<int> genTtbarIdHandle;
    iEvent.getByToken(genTtbarIdToken_, genTtbarIdHandle);
    b_genTtbarId = *genTtbarIdHandle;

    edm::Handle<int> NgenJet30Handle;
    iEvent.getByToken(NgenJet30Token_, NgenJet30Handle);
    b_NgenJet30 = *NgenJet30Handle;

    edm::Handle<int> genTtbarLeptonDecayHandle;
    iEvent.getByToken(genTtbarLeptonDecayToken_, genTtbarLeptonDecayHandle);
    b_genTtbarLeptonDecay = *genTtbarLeptonDecayHandle;

  }
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    ttree_->Fill();
    return;
  }
  cutflow_[1][b_channel]++;

  // const reco::Vertex &PV = vertices->front();
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  b_nvertex = *nGoodVertexHandle;

  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;
  // if (!b_filtered){
  //   ttree_->Fill();
  //   return;
  // }
  cutflow_[2][b_channel]++;
  
  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);
  
  // Find leptons and sort by pT
  LeptonCollection recolep;
  selectMuons(*muons, recolep);
  selectElecs(*electrons, recolep);
  if (recolep.size() < 2){
    ttree_->Fill();
    return;
  }
  cutflow_[3][b_channel]++;

  sort(recolep.begin(), recolep.end(), GtByCandPt());
  int lep1_idx=0, lep2_idx=1;
  //for emu
  if(std::abs(recolep[1].pdgId())==11 && std::abs(recolep[0].pdgId())==13){
     lep1_idx = 1; lep2_idx = 0;
  }

  const cat::Lepton& recolep1 = recolep[lep1_idx];
  const cat::Lepton& recolep2 = recolep[lep2_idx];

  // Determine channel
  const int pdgIdSum = std::abs(recolep1.pdgId()) + std::abs(recolep2.pdgId());
  if (pdgIdSum == 24) b_channel = CH_MUEL; // emu
  if (pdgIdSum == 22) b_channel = CH_ELEL; // ee
  if (pdgIdSum == 26) b_channel = CH_MUMU; // mumu

  // Trigger results
  edm::Handle<int> trigHandle;
  if      ( b_channel == CH_ELEL ) iEvent.getByToken(trigTokenELEL_, trigHandle);
  else if ( b_channel == CH_MUMU ) iEvent.getByToken(trigTokenMUMU_, trigHandle);
  else if ( b_channel == CH_MUEL ) iEvent.getByToken(trigTokenMUEL_, trigHandle);
  b_tri = *trigHandle;

  b_lep1_pt = recolep1.pt(); b_lep1_eta = recolep1.eta(); b_lep1_phi = recolep1.phi();
  b_lep2_pt = recolep2.pt(); b_lep2_eta = recolep2.eta(); b_lep2_phi = recolep2.phi();

  if (pdgIdSum == 24) {  b_lep1_RelIso = recolep1.relIso();     b_lep2_RelIso = recolep2.relIso(0.4); } // emu
  if (pdgIdSum == 22) {  b_lep1_RelIso = recolep1.relIso();     b_lep2_RelIso = recolep2.relIso();    } // ee
  if (pdgIdSum == 26) {  b_lep1_RelIso = recolep1.relIso(0.4);  b_lep2_RelIso = recolep2.relIso(0.4); } // mumu

  const auto tlv_ll = recolep1.p4()+recolep2.p4();
  b_ll_pt = tlv_ll.Pt(); b_ll_eta = tlv_ll.Eta(); b_ll_phi = tlv_ll.Phi(); b_ll_m = tlv_ll.M();

  if (b_ll_m < 20. || recolep1.charge() * recolep2.charge() > 0){
    ttree_->Fill();
    return;
  }
  else b_step1 = true;
  b_step = 1;
  cutflow_[4][b_channel]++;

  if ( (b_channel == CH_MUEL) || ((b_ll_m < 76) || (b_ll_m > 106)) ){
    b_step2 = true;
    b_step = 2;
    cutflow_[5][b_channel]++;
  }

  JetCollection&& selectedJets = selectJets(*jets, recolep);
  JetCollection&& selectedBJetsL = selectBJets(selectedJets,0.605);
  JetCollection&& selectedBJetsM = selectBJets(selectedJets,0.890);
  JetCollection&& selectedBJetsT = selectBJets(selectedJets,0.970);

  int idx=0;
  std::map<int,double> mapJetBDiscriminator;
  for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
    double bDisCSV= jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    int flavor = jet1->partonFlavour();
    mapJetBDiscriminator[idx] = bDisCSV;
    idx++;
    b_jets_pt.push_back(jit->p4().pt());
    b_jets_eta.push_back(jit->p4().eta());
    b_jets_phi.push_back(jit->p4().phi());
    b_jets_flavor.push_back(flavor);
    b_jets_bDiscriminatorCSV.push_back(bDiscriminator);
  }
  if (runOnMC_){
     double csvWgtHF, csvWgtLF, csvWgtCF;
     b_csvweight         = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF); 
     b_csvweight_JES_Up  = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_JES_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_LF_Up   = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_LF_Down = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_HF_Up   = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_HF_Down = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_HF_Stats1_Up  = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_HF_Stats1_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_HF_Stats2_Up  = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_HF_Stats2_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_LF_Stats1_Up  = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_LF_Stats1_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_LF_Stats2_Up  = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_LF_Stats2_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_Charm_Err1_Up = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_Charm_Err1_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_Charm_Err2_Up  = csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
     b_csvweight_Charm_Err2_Down= csvWeight->getCSVWeight(b_jets_pt, b_jets_eta, b_jets_phi, b_jets_flavor, 0, csvWgtHF, csvWgtLF, csvWgtCF);
  }
  //csvd order
  std::vector< std::pair<int,double> > vecJetBDisc(mapJetBDiscriminator.begin(), mapJetBDiscriminator.end());
  std::sort(vecJetBDisc.begin(), vecJetBDisc.end(), bigger_second<data_t>());
  for(std::vector< std::pair<int,double> >::iterator it = vecJetBDisc.begin() ; it != vecJetBDisc.end(); ++it)
      b_csvd_jetid.push_back((*it).first)

  const auto met = mets->front().p4();
  b_met = met.pt();
  b_metphi = met.phi();
  b_njet30 = selectedJets.size();
  b_nbjetM30 = selectedBJetsM.size();
  b_nbjetMT = selectedBJetsT.size();

  if (selectedJets.size() >1 ){
    b_step3 = true;
    if (b_step == 2){
      ++b_step;
      cutflow_[6][b_channel]++;
    }
  }

  if ((b_channel == CH_MUEL) || (b_met > 40.)){
    b_step4 = true;
    if (b_step == 3){
      ++b_step;
      cutflow_[7][b_channel]++;
    }
  }
  
  if (selectedBJets.size() > 0){
    b_step5 = true;
    if (b_step == 4){
      ++b_step;
      cutflow_[8][b_channel]++;
    }
  }

/*  ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
  //int kin=0;
  math::XYZTLorentzVector top1, top2, nu1, nu2;
  double maxweight=0;
  //const cat::Jet* kinj1, * kinj2;

  const auto recolepLV1= recolep1.p4();
  const auto recolepLV2= recolep2.p4();
  math::XYZTLorentzVector inputLV[5] = {met, recolepLV1, recolepLV2};

  for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
    const auto recojet1= jet1->p4();
    for (auto jet2 = next(jet1); jet2 != end; ++jet2){

      const auto recojet2= jet2->p4();

      b_jet1_pt = recojet1.Pt();
      b_jet1_eta = recojet1.Eta();
      b_jet2_pt = recojet2.Pt();
      b_jet2_eta = recojet2.Eta();

      inputLV[3] = recojet1;
      inputLV[4] = recojet2;
      solver_->solve(inputLV);
      const double weight1 = solver_->quality();
      inputLV[3] = recojet2;
      inputLV[4] = recojet1;
      solver_->solve(inputLV);
      const double weight2 = solver_->quality();

      if ( weight2 > maxweight and weight2 >= weight1 ) {
        nu1 = solver_->nu1();
        nu2 = solver_->nu2();
        maxweight = weight2;
      }
      else if ( weight1 > maxweight and weight1 >= weight2 ) {
        // Re-solve with previous jet combinations
        // Weights are re-calculated since there can be very little difference due to random number effect in smearing algorithm
	inputLV[3] = recojet1;
        inputLV[4] = recojet2;
        solver_->solve(inputLV);
        nu1 = solver_->nu1();
        nu2 = solver_->nu2();
        maxweight = solver_->quality();
      }
      else continue;

      top1 = recolepLV1+recojet1+nu1;
      top2 = recolepLV2+recojet2+nu2;
    }
  }

  b_top1_pt = top1.Pt();
  b_top1_eta = top1.Eta();
  b_top1_phi = top1.Phi();
  b_top1_rapi = top1.Rapidity();
  b_top2_pt = top2.Pt();
  b_top2_eta = top2.Eta();
  b_top2_phi = top2.Phi();
  b_top2_rapi = top2.Rapidity();

  auto ttbar = top1+top2;
  b_ttbar_pt = ttbar.Pt();
  b_ttbar_eta = ttbar.Eta();
  b_ttbar_phi = ttbar.Phi();
  b_ttbar_m = ttbar.M();
  b_ttbar_rapi = ttbar.Rapidity();

  b_maxweight = maxweight;
  if (maxweight){
    b_step6 = true;
    if (b_step == 5){
      ++b_step;
    }
  }
*/
  //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
  // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_met, b_ll_mass, b_maxweight);
  ttree_->Fill();
}

const reco::Candidate* TtbarBbbarDiLeptonAnalyzer::getLast(const reco::Candidate* p) const
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
    {
      const reco::Candidate* dau = p->daughter(i);
      if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
    }
  return p;
}

void TtbarBbbarDiLeptonAnalyzer::selectMuons(const cat::MuonCollection& muons, LeptonCollection& selmuons) const
{
  for (auto& mu : muons) {
    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    //if (!mu.isMediumMuon()) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.12) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }
}

void TtbarBbbarDiLeptonAnalyzer::selectElecs(const cat::ElectronCollection& elecs, LeptonCollection& selelecs) const
{
  for (auto& el : elecs) {
    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    //if (!el.electronID("cutBasedElectronID-Spring15-50ns-V1-standalone-medium")) continue;
    //if (el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium") == 0) continue;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) continue;
    //if (!el.passConversionVeto()) continue;
    //if (!el.isPF()) continue;

    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
}

cat::JetCollection TtbarBbbarDiLeptonAnalyzer::selectJets(const cat::JetCollection& jets, const LeptonCollection& recolep) const
{
  cat::JetCollection seljets;
  for (auto& jet : jets) {
    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 2.4)  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep.p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}

cat::JetCollection TtbarBbbarDiLeptonAnalyzer::selectBJets(const JetCollection& jets, double workingpoint) const
{
  //https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV 
  // 25ns pfCombinedInclusiveSecondaryVertexV2BJetTags :L,M,T 0.605, 0.890, 0.970
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < workingpoint) continue;
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarBbbarDiLeptonAnalyzer);
