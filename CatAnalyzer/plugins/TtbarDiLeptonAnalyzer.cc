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

#include "CATTools/CommonTools/interface/ScaleFactorGetter.h"
//#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TTree.h"
//#include "TLorentzVector.h"

using namespace std;
using namespace cat;

class TtbarDiLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer();

  enum {
    CH_NONE=0, CH_MUEL=1, CH_ELEL=2, CH_MUMU=3
  };
  enum sys_e {sys_nom, 
    sys_jes_u, sys_jes_d, sys_jer_u, sys_jer_d,
    sys_mu_u, sys_mu_d, sys_el_u, sys_el_d,
    sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
    sys_btag_u, sys_btag_d,
    nsys_e
  };

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  void resetBr();
  float selectMuons(const cat::MuonCollection& muons, ParticleCollection& selmuons, sys_e sys) const;
  float selectElecs(const cat::ElectronCollection& elecs, ParticleCollection& selelecs, sys_e sys) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const ParticleCollection& recolep, sys_e sys);
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  const reco::Candidate* getLast(const reco::Candidate* p) const;

  ScaleFactorGetter muonSF_, elecSF_;
  float getSF(const cat::Particle& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 13 ) {
      const double pt = p.pt(), eta = std::abs(p.eta());
      if      ( sys == sys_mueff_u ) return muonSF_(eta, pt,  1);
      else if ( sys == sys_mueff_d ) return muonSF_(eta, pt, -1);
      else return muonSF_(eta, pt, 0);
    }
    else {
      const double pt = p.pt(), eta = p.eta();
      if      ( sys == sys_eleff_u ) return elecSF_(eta, pt,  1);
      else if ( sys == sys_eleff_d ) return elecSF_(eta, pt, -1);
      else return elecSF_(eta, pt, 0);
    }
    return 1;
  }

  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genweightToken_, puweightToken_, puweightToken_up_, puweightToken_dn_;
  edm::EDGetTokenT<vector<float>> pdfweightToken_;
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<vector<int> > partonTop_modes_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTop_genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > pseudoTop_leptons_, pseudoTop_neutrinos_, pseudoTop_jets_;

  std::vector<TTree*> ttree_;
  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_tri, b_filtered;
  float b_met, b_weight, b_puweight, b_puweight_up, b_puweight_dn, b_genweight, b_lepweight, b_btagweight, b_btagweight_up, b_btagweight_dn;
  std::vector<float> b_pdfWeights;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;

  float b_partontop1_pt, b_partontop1_eta, b_partontop1_phi, b_partontop1_rapi, b_partontop1_m;
  float b_partontop2_pt, b_partontop2_eta, b_partontop2_phi, b_partontop2_rapi, b_partontop2_m;
  float b_partonttbar_pt, b_partonttbar_eta, b_partonttbar_dphi, b_partonttbar_m, b_partonttbar_rapi;

  int b_partonChannel, b_partonMode1, b_partonMode2;
  float b_partonlep1_pt, b_partonlep1_eta;
  float b_partonlep2_pt, b_partonlep2_eta;
  float b_partonjet1_pt, b_partonjet1_eta;
  float b_partonjet2_pt, b_partonjet2_eta;
  bool b_partonInPhase, b_partonInPhaseJet, b_partonInPhaseLep;

  float b_gentop1_pt, b_gentop1_eta, b_gentop1_phi, b_gentop1_rapi, b_gentop1_m;
  float b_gentop2_pt, b_gentop2_eta, b_gentop2_phi, b_gentop2_rapi, b_gentop2_m;
  float b_genttbar_pt, b_genttbar_eta, b_genttbar_dphi, b_genttbar_m, b_genttbar_rapi;

  int b_pseudoTopChannel;
  float b_genlep1_pt, b_genlep1_eta;
  float b_genlep2_pt, b_genlep2_eta;
  float b_genjet1_pt, b_genjet1_eta;
  float b_genjet2_pt, b_genjet2_eta;
  bool b_pseudoInPhase;

  float b_jet1_pt, b_jet1_eta, b_jet1_CSVInclV2;
  float b_jet2_pt, b_jet2_eta, b_jet2_CSVInclV2;
  float b_top1_pt, b_top1_eta, b_top1_phi, b_top1_rapi, b_top1_m;
  float b_top2_pt, b_top2_eta, b_top2_phi, b_top2_rapi, b_top2_m;
  float b_ttbar_pt, b_ttbar_eta, b_ttbar_dphi, b_ttbar_m, b_ttbar_rapi;
  float b_maxweight;
  int b_is3lep;

  //std::unique_ptr<TtFullLepKinSolver> solver;
  std::unique_ptr<KinematicSolver> solver_;
  //enum TTbarMode { CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON };
  //enum DecayMode { CH_HADRON = 1, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

  const static int NCutflow = 10;
  std::vector<std::vector<int> > cutflow_;
};
//
// constructors and destructor
//
TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  pdfweightToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("pdfweight"));
  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightToken_up_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_up"));
  puweightToken_dn_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_dn"));
  trigTokenMUEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));

  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  typedef std::vector<double> vdouble;

  const auto muonSet = iConfig.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("effSF");
  muonSF_.set(muonSFSet.getParameter<vdouble>("abseta_bins"),
              muonSFSet.getParameter<vdouble>("pt_bins"),
              muonSFSet.getParameter<vdouble>("values"),
              muonSFSet.getParameter<vdouble>("errors"));

  const auto elecSet = iConfig.getParameter<edm::ParameterSet>("electron");
  elecToken_ = consumes<cat::ElectronCollection>(elecSet.getParameter<edm::InputTag>("src"));
  const auto elecSFSet = elecSet.getParameter<edm::ParameterSet>("effSF");
  elecSF_.set(elecSFSet.getParameter<vdouble>("pt_bins"),
              elecSFSet.getParameter<vdouble>("eta_bins"),
              elecSFSet.getParameter<vdouble>("values"),
              elecSFSet.getParameter<vdouble>("errors"));

  partonTop_channel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("partonTop_channel"));
  partonTop_modes_   = consumes<vector<int> >(iConfig.getParameter<edm::InputTag>("partonTop_modes"));
  partonTop_genParticles_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("partonTop_genParticles"));
  pseudoTop_leptons_   = consumes<edm::View<reco::Candidate> >(edm::InputTag("pseudoTop", "leptons"));
  pseudoTop_neutrinos_   = consumes<edm::View<reco::Candidate> >(edm::InputTag("pseudoTop", "neutrinos"));
  pseudoTop_jets_ = consumes<edm::View<reco::Candidate> >(edm::InputTag("pseudoTop", "jets"));

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

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  const std::string sys_name[nsys_e] = {
    "nom",
    "jes_u", "jes_d", "jer_u", "jer_d",
    "mu_u", "mu_d", "el_u", "el_d",
    "mueff_u", "mueff_d", "eleff_u", "eleff_d",
    "btag_u", "btag_d"
  };
  for (int sys = 0; sys < nsys_e; ++sys){
    ttree_.push_back(fs->make<TTree>(sys_name[sys].c_str(), sys_name[sys].c_str()));
    auto tr = ttree_.back();
    tr->Branch("nvertex", &b_nvertex, "nvertex/I");
    tr->Branch("step", &b_step, "step/I");
    tr->Branch("channel", &b_channel, "channel/I");
    tr->Branch("njet", &b_njet, "njet/I");
    tr->Branch("nbjet", &b_nbjet, "nbjet/I");
    tr->Branch("step1", &b_step1, "step1/O");
    tr->Branch("step2", &b_step2, "step2/O");
    tr->Branch("step3", &b_step3, "step3/O");
    tr->Branch("step4", &b_step4, "step4/O");
    tr->Branch("step5", &b_step5, "step5/O");
    tr->Branch("step6", &b_step6, "step6/O");
    tr->Branch("tri", &b_tri, "tri/O");
    tr->Branch("filtered", &b_filtered, "filtered/O");
    tr->Branch("met", &b_met, "met/F");
    tr->Branch("weight", &b_weight, "weight/F");
    tr->Branch("pdfWeihgts","std::vector<float>",&b_pdfWeights);
    tr->Branch("puweight", &b_puweight, "puweight/F");
    tr->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
    tr->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
    tr->Branch("genweight", &b_genweight, "genweight/F");
    tr->Branch("lepweight", &b_lepweight, "lepweight/F");
    tr->Branch("btagweight", &b_btagweight, "btagweight/F");
    tr->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
    tr->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");

    tr->Branch("lep1_pt", &b_lep1_pt, "lep1_pt/F");
    tr->Branch("lep1_eta", &b_lep1_eta, "lep1_eta/F");
    tr->Branch("lep1_phi", &b_lep1_phi, "lep1_phi/F");
    tr->Branch("lep2_pt", &b_lep2_pt, "lep2_pt/F");
    tr->Branch("lep2_eta", &b_lep2_eta, "lep2_eta/F");
    tr->Branch("lep2_phi", &b_lep2_phi, "lep2_phi/F");
    tr->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
    tr->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
    tr->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
    tr->Branch("ll_m", &b_ll_m, "ll_m/F");

    tr->Branch("partontop1_pt", &b_partontop1_pt, "partontop1_pt/F");
    tr->Branch("partontop1_eta", &b_partontop1_eta, "partontop1_eta/F");
    tr->Branch("partontop1_rapi", &b_partontop1_rapi, "partontop1_rapi/F");
    tr->Branch("partontop1_m", &b_partontop1_m, "partontop1_m/F");
    tr->Branch("partontop2_pt", &b_partontop2_pt, "partontop2_pt/F");
    tr->Branch("partontop2_eta", &b_partontop2_eta, "partontop2_eta/F");
    tr->Branch("partontop2_rapi", &b_partontop2_rapi, "partontop2_rapi/F");
    tr->Branch("partontop2_m", &b_partontop2_m, "partontop2_m/F");
    tr->Branch("partonttbar_pt", &b_partonttbar_pt, "partonttbar_pt/F");
    tr->Branch("partonttbar_eta", &b_partonttbar_eta, "partonttbar_eta/F");
    tr->Branch("partonttbar_dphi", &b_partonttbar_dphi, "partonttbar_dphi/F");
    tr->Branch("partonttbar_rapi", &b_partonttbar_rapi, "partonttbar_rapi/F");
    tr->Branch("partonttbar_m", &b_partonttbar_m, "partonttbar_m/F");

    tr->Branch("parton_channel", &b_partonChannel, "parton_channel/I");
    tr->Branch("parton_mode1", &b_partonMode1, "parton_mode1/I");
    tr->Branch("partonlep1_pt", &b_partonlep1_pt, "partonlep1_pt/F");
    tr->Branch("partonlep1_eta", &b_partonlep1_eta, "partonlep1_eta/F");
    tr->Branch("partonlep2_pt", &b_partonlep2_pt, "partonlep2_pt/F");
    tr->Branch("partonlep2_eta", &b_partonlep2_eta, "partonlep2_eta/F");
    tr->Branch("partonjet1_pt", &b_partonjet1_pt, "partonjet1_pt/F");
    tr->Branch("partonjet1_eta", &b_partonjet1_eta, "partonjet1_eta/F");
    tr->Branch("partonjet2_pt", &b_partonjet2_pt, "partonjet2_pt/F");
    tr->Branch("partonjet2_eta", &b_partonjet2_eta, "partonjet2_eta/F");
    tr->Branch("parton_mode2", &b_partonMode2, "parton_mode2/I");
    tr->Branch("partonInPhase", &b_partonInPhase, "partonInPhase/O");
    tr->Branch("partonInPhaseLep", &b_partonInPhaseLep, "partonInPhaseLep/O");
    tr->Branch("partonInPhaseJet", &b_partonInPhaseJet, "partonInPhaseJet/O");

    tr->Branch("gentop1_pt", &b_gentop1_pt, "gentop1_pt/F");
    tr->Branch("gentop1_eta", &b_gentop1_eta, "gentop1_eta/F");
    tr->Branch("gentop1_rapi", &b_gentop1_rapi, "gentop1_rapi/F");
    tr->Branch("gentop1_m", &b_gentop1_m, "gentop1_m/F");
    tr->Branch("gentop2_pt", &b_gentop2_pt, "gentop2_pt/F");
    tr->Branch("gentop2_eta", &b_gentop2_eta, "gentop2_eta/F");
    tr->Branch("gentop2_rapi", &b_gentop2_rapi, "gentop2_rapi/F");
    tr->Branch("gentop2_m", &b_gentop2_m, "gentop2_m/F");
    tr->Branch("genttbar_pt", &b_genttbar_pt, "genttbar_pt/F");
    tr->Branch("genttbar_eta", &b_genttbar_eta, "genttbar_eta/F");
    tr->Branch("genttbar_dphi", &b_genttbar_dphi, "genttbar_dphi/F");
    tr->Branch("genttbar_rapi", &b_genttbar_rapi, "genttbar_rapi/F");
    tr->Branch("genttbar_m", &b_genttbar_m, "genttbar_m/F");

    tr->Branch("pseudoTop_channel", &b_pseudoTopChannel, "pseudoTop_channel/I");
    tr->Branch("genlep1_pt", &b_genlep1_pt, "genlep1_pt/F");
    tr->Branch("genlep1_eta", &b_genlep1_eta, "genlep1_eta/F");
    tr->Branch("genlep2_pt", &b_genlep2_pt, "genlep2_pt/F");
    tr->Branch("genlep2_eta", &b_genlep2_eta, "genlep2_eta/F");
    tr->Branch("genjet1_pt", &b_genjet1_pt, "genjet1_pt/F");
    tr->Branch("genjet1_eta", &b_genjet1_eta, "genjet1_eta/F");
    tr->Branch("genjet2_pt", &b_genjet2_pt, "genjet2_pt/F");
    tr->Branch("genjet2_eta", &b_genjet2_eta, "genjet2_eta/F");
    tr->Branch("pseudoInPhase", &b_pseudoInPhase, "pseudoInPhase/O");

    tr->Branch("jet1_pt", &b_jet1_pt, "jet1_pt/F");
    tr->Branch("jet2_pt", &b_jet2_pt, "jet2_pt/F");
    tr->Branch("jet1_eta", &b_jet1_eta, "jet1_eta/F");
    tr->Branch("jet2_eta", &b_jet2_eta, "jet2_eta/F");
    tr->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
    tr->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");

    tr->Branch("top1_pt", &b_top1_pt, "top1_pt/F");
    tr->Branch("top1_eta", &b_top1_eta, "top1_eta/F");
    tr->Branch("top1_phi", &b_top1_phi, "top1_phi/F");
    tr->Branch("top1_rapi", &b_top1_rapi, "top1_rapi/F");
    tr->Branch("top1_m", &b_top1_m, "top1_m/F");
    tr->Branch("top2_pt", &b_top2_pt, "top2_pt/F");
    tr->Branch("top2_eta", &b_top2_eta, "top2_eta/F");
    tr->Branch("top2_phi", &b_top2_phi, "top2_phi/F");
    tr->Branch("top2_rapi", &b_top2_rapi, "top2_rapi/F");
    tr->Branch("top2_m", &b_top2_m, "top2_m/F");
    tr->Branch("ttbar_pt", &b_ttbar_pt, "ttbar_pt/F");
    tr->Branch("ttbar_eta", &b_ttbar_eta, "ttbar_eta/F");
    tr->Branch("ttbar_dphi", &b_ttbar_dphi, "ttbar_dphi/F");
    tr->Branch("ttbar_rapi", &b_ttbar_rapi, "ttbar_rapi/F");
    tr->Branch("ttbar_m", &b_ttbar_m, "ttbar_m/F");

    tr->Branch("is3lep", &b_is3lep, "is3lep/I");
  }

  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});
}

TtbarDiLeptonAnalyzer::~TtbarDiLeptonAnalyzer()
{
  cout <<"     cut flow   emu    ee    mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step "<< i << "    "<< cutflow_[i][0] <<  "   "<< cutflow_[i][1] << "   " << cutflow_[i][2] << "   " << cutflow_[i][3]<< endl;
  }
}

void TtbarDiLeptonAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&)
{
  if ( dynamic_cast<DESYSmearedSolver*>(solver_.get()) != 0 ) {
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(lumi.index());
    dynamic_cast<DESYSmearedSolver*>(solver_.get())->setRandom(&engine);
  }
}

void TtbarDiLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const bool runOnMC = !iEvent.isRealData();

  cutflow_[0][0]++;

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();

    edm::Handle<int> partonTop_channel;
    if ( iEvent.getByToken(partonTop_channel_, partonTop_channel)){
      
      if (sys == sys_nom){
	edm::Handle<vector<float>> pdfweightHandle;
	iEvent.getByToken(pdfweightToken_, pdfweightHandle);
	for (const float & aPdfWeight : *pdfweightHandle){
	  b_pdfWeights.push_back(aPdfWeight);  
	}
      }
      
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
        b_partontop1_pt = parton1->pt();
        b_partontop1_eta = parton1->eta();
        b_partontop1_phi = parton1->phi();
        b_partontop1_rapi = parton1->rapidity();
        b_partontop1_m = parton1->mass();
        b_partontop2_pt = parton2->pt();
        b_partontop2_eta = parton2->eta();
        b_partontop2_rapi = parton2->rapidity();
        b_partontop2_m = parton2->mass();

        // Get TTbar
        auto partonttbar = parton1->p4()+parton2->p4();
        b_partonttbar_pt = partonttbar.Pt();
        b_partonttbar_eta = partonttbar.Eta();
        b_partonttbar_m = partonttbar.M();
        b_partonttbar_rapi = partonttbar.Rapidity();
        b_partonttbar_dphi = deltaPhi((parton1->p4()).Phi(), (parton2->p4()).Phi());

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
            // Fill lepton informations
            b_partonjet1_pt = partonB1->pt();
            b_partonjet1_eta = partonB1->eta();
            b_partonjet2_pt = partonB2->pt();
            b_partonjet2_eta = partonB2->eta();
          }
        }
        if (b_partonInPhaseJet && b_partonInPhaseLep) b_partonInPhase = true;
      }

      // Start to build pseudo top
      b_pseudoTopChannel = CH_NONE;
      edm::Handle<edm::View<reco::Candidate> > pseudoTopLeptonHandle;
      edm::Handle<edm::View<reco::Candidate> > pseudoTopNeutrinoHandle;
      edm::Handle<edm::View<reco::Candidate> > pseudoTopJetHandle;
      iEvent.getByToken(pseudoTop_leptons_, pseudoTopLeptonHandle);
      iEvent.getByToken(pseudoTop_neutrinos_, pseudoTopNeutrinoHandle);
      iEvent.getByToken(pseudoTop_jets_, pseudoTopJetHandle);
      do {
        // Basic lepton, jet multiplicity
        if ( pseudoTopLeptonHandle->size() < 2 or pseudoTopJetHandle->size() < 2 or pseudoTopNeutrinoHandle->size() < 2 ) break;

        std::vector<size_t> leptonIdxs, neutrinoIdxs, bjetIdxs;
        // Lepton acceptance cuts
        for ( size_t i=0, n=pseudoTopLeptonHandle->size(); i<n; ++i ) {
          const auto& x = pseudoTopLeptonHandle->at(i);
          if ( x.pt() < 20 or std::abs(x.eta()) > 2.4 ) continue;
          if ( abs(x.pdgId()) != 11 and abs(x.pdgId()) != 13 ) continue;
          leptonIdxs.push_back(i);
        }
        if ( leptonIdxs.size() < 2 ) break;
        std::nth_element(leptonIdxs.begin(), leptonIdxs.begin()+2, leptonIdxs.end(),
            [&](size_t i, size_t j){return pseudoTopLeptonHandle->at(i).pt() > pseudoTopLeptonHandle->at(j).pt();});
        const auto lepton1 = pseudoTopLeptonHandle->at(leptonIdxs[0]).p4();
        const auto lepton2 = pseudoTopLeptonHandle->at(leptonIdxs[1]).p4();
        const int pseudoW1DauId = abs(pseudoTopLeptonHandle->at(leptonIdxs[0]).pdgId());
        const int pseudoW2DauId = abs(pseudoTopLeptonHandle->at(leptonIdxs[1]).pdgId());
        switch ( pseudoW1DauId+pseudoW2DauId ) {
          case 22: b_pseudoTopChannel = CH_ELEL; break;
          case 26: b_pseudoTopChannel = CH_MUMU; break;
          case 24: b_pseudoTopChannel = CH_MUEL; break;
          default: b_pseudoTopChannel = CH_NONE;
        }

        //std::nth_element(neutrinoIdxs.begin(), neutrinoIdxs.begin()+2, neutrinoIdxs.end(),
        //                 [&](size_t i, size_t j){return pseudoTopLeptonHandle->at(i).pt() > pseudoTopLeptonHandle->at(j).pt();});
        auto nu1 = pseudoTopNeutrinoHandle->at(0).p4(), nu2 = pseudoTopNeutrinoHandle->at(1).p4();

        // Jet acceptance and generator level b tag
        for ( size_t i=0, n=pseudoTopJetHandle->size(); i<n; ++i ) {
          const auto& x = pseudoTopJetHandle->at(i);
          if ( x.pt() < 30 or std::abs(x.eta()) > 2.4 ) continue;
          if ( abs(x.pdgId()) != 5 ) continue;
          bjetIdxs.push_back(i);
        }
        if ( bjetIdxs.size() < 2 ) break;
        std::nth_element(bjetIdxs.begin(), bjetIdxs.begin()+2, bjetIdxs.end(),
            [&](size_t i, size_t j){return pseudoTopJetHandle->at(i).pt() > pseudoTopJetHandle->at(j).pt();});
        auto bjet1 = pseudoTopJetHandle->at(bjetIdxs[0]).p4(), bjet2 = pseudoTopJetHandle->at(bjetIdxs[1]).p4();

        // Do the W combinations
        auto w1 = lepton1 + nu1;
        auto w2 = lepton2 + nu2;
        if ( true ) {
          const auto w1Alt = lepton1 + nu2;
          const auto w2Alt = lepton2 + nu1;

          const double wMass = 80.4;
          const double dm = std::abs(w1.mass()-wMass)+std::abs(w2.mass()-wMass);
          const double dmAlt = std::abs(w1Alt.mass()-wMass)+std::abs(w2Alt.mass()-wMass);
          if ( dm > dmAlt ) { w1 = w1Alt; w2 = w2Alt; std::swap(nu1, nu2); }
        }
        // Do the top combinations
        auto gentop1 = w1 + bjet1;
        auto gentop2 = w2 + bjet2;
        if ( true ) {
          const auto t1Alt = w1 + bjet2;
          const auto t2Alt = w2 + bjet1;

          const double tMass = 172.5;
          const double dm = std::abs(gentop1.mass()-tMass)+std::abs(gentop2.mass()-tMass);
          const double dmAlt = std::abs(t1Alt.mass()-tMass)+std::abs(t2Alt.mass()-tMass);
          if ( dm > dmAlt ) { gentop1 = t1Alt; gentop2 = t2Alt; std::swap(bjet1, bjet2); }
        }

        if (gentop1.Pt() < gentop2.Pt()) { swap(gentop1, gentop2); }
        b_gentop1_pt = gentop1.Pt();
        b_gentop1_eta = gentop1.Eta();
        b_gentop1_phi = gentop1.Phi();
        b_gentop1_rapi = gentop1.Rapidity();
        b_gentop1_m = gentop1.M();
        b_gentop2_pt = gentop2.Pt();
        b_gentop2_eta = gentop2.Eta();
        b_gentop2_phi = gentop2.Phi();
        b_gentop2_rapi = gentop2.Rapidity();
        b_gentop2_m = gentop2.M();

        // Get Top quark pairs
        auto genttbar = gentop1+gentop2;
        b_genttbar_pt = genttbar.Pt();
        b_genttbar_eta = genttbar.Eta();
        b_genttbar_dphi = deltaPhi(gentop1.Phi(), gentop2.Phi());
        b_genttbar_m = genttbar.M();
        b_genttbar_rapi = genttbar.Rapidity();

        b_genlep1_pt = lepton1.pt();
        b_genlep1_eta = lepton1.eta();
        b_genlep2_pt = lepton2.pt();
        b_genlep2_eta = lepton2.eta();

        if (bjet1.Pt() < bjet2.Pt()) { swap(bjet1, bjet2); }
        b_genjet1_pt = bjet1.pt();
        b_genjet1_eta = bjet1.eta();
        b_genjet2_pt = bjet2.pt();
        b_genjet2_eta = bjet2.eta();

      } while ( false );
    }
    //  } //sys == sys_nom
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
      iEvent.getByToken(genweightToken_, genweightHandle);
      b_genweight = (*genweightHandle);
      b_weight = b_genweight*b_puweight;
    }

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()){ // skip the event if no PV found
      ttree_[sys]->Fill();
      continue;
    }
    cutflow_[1][b_channel]++;

    // const reco::Vertex &PV = vertices->front();
    edm::Handle<int> nGoodVertexHandle;
    iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
    b_nvertex = *nGoodVertexHandle;

    edm::Handle<int> lumiSelectionHandle;
    iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
    if (!runOnMC){
      if (*lumiSelectionHandle == 0) return;
    }

    edm::Handle<int> recoFiltersHandle;
    iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
    b_filtered = *recoFiltersHandle == 0 ? false : true;
    // if (!b_filtered){
    //   ttree_[sys]->Fill();
    //   continue;
    // }
    cutflow_[2][b_channel]++;

    edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
    edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
    edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
    edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);

    // Find leptons and sort by pT
    ParticleCollection recolep;
    selectMuons(*muons, recolep, (sys_e)sys);
    selectElecs(*electrons, recolep, (sys_e)sys);
    //b_lepweight = mu_weight * el_weight;
    if (recolep.size() < 2){
      ttree_[sys]->Fill();
      continue;
    }
    cutflow_[3][b_channel]++;

    sort(recolep.begin(), recolep.end(), GtByCandPt());
    recolep.erase(recolep.begin()+2,recolep.end());
    const cat::Particle& recolep1 = recolep[0];
    const cat::Particle& recolep2 = recolep[1];

    // Determine channel
    const int pdgIdSum = std::abs(recolep1.pdgId()) + std::abs(recolep2.pdgId());
    if (pdgIdSum == 24) b_channel = CH_MUEL; // emu
    if (pdgIdSum == 22) b_channel = CH_ELEL; // ee
    if (pdgIdSum == 26) b_channel = CH_MUMU; // mumu

    b_lepweight = getSF(recolep1, sys)*getSF(recolep2, sys);
    b_weight *= b_lepweight;

    // Trigger results
    edm::Handle<int> trigHandle;
    if      ( b_channel == CH_ELEL ) iEvent.getByToken(trigTokenELEL_, trigHandle);
    else if ( b_channel == CH_MUMU ) iEvent.getByToken(trigTokenMUMU_, trigHandle);
    else if ( b_channel == CH_MUEL ) iEvent.getByToken(trigTokenMUEL_, trigHandle);
    b_tri = *trigHandle;

    b_lep1_pt = recolep1.pt(); b_lep1_eta = recolep1.eta(); b_lep1_phi = recolep1.phi();
    b_lep2_pt = recolep2.pt(); b_lep2_eta = recolep2.eta(); b_lep2_phi = recolep2.phi();
    const auto tlv_ll = recolep1.p4()+recolep2.p4();
    b_ll_pt = tlv_ll.Pt(); b_ll_eta = tlv_ll.Eta(); b_ll_phi = tlv_ll.Phi(); b_ll_m = tlv_ll.M();

    if (b_ll_m < 20. || recolep1.charge() * recolep2.charge() > 0){
      ttree_[sys]->Fill();
      continue;
    }
    b_step1 = true;
    b_step = 1;
    cutflow_[4][b_channel]++;

    if ( (b_channel == CH_MUEL) || ((b_ll_m < 76) || (b_ll_m > 106)) ){
      b_step2 = true;
      b_step = 2;
      cutflow_[5][b_channel]++;
    }

    JetCollection&& selectedJets = selectJets(*jets, recolep, (sys_e)sys);
    JetCollection&& selectedBJets = selectBJets(selectedJets);

    const auto met = mets->front().p4();
    b_met = met.pt();
    b_njet = selectedJets.size();
    b_nbjet = selectedBJets.size();

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

    ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
    //int kin=0;
    math::XYZTLorentzVector top1, top2, nu1, nu2, bjet1, bjet2;
    double maxweight=0;
    //const cat::Jet* kinj1, * kinj2;

    const auto recolepLV1= recolep1.p4();
    const auto recolepLV2= recolep2.p4();
    math::XYZTLorentzVector inputLV[5] = {met, recolepLV1, recolepLV2};

    for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
      const auto recojet1= jet1->p4();
      for (auto jet2 = next(jet1); jet2 != end; ++jet2){

        const auto recojet2= jet2->p4();

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

        //saving results
        bjet1 = recojet1;
        bjet2 = recojet2;

        top1 = recolepLV1+inputLV[3]+nu1;
        top2 = recolepLV2+inputLV[4]+nu2;

      }
    }

    if (bjet1.Pt() < bjet2.Pt()) { swap(bjet1, bjet2); }
    b_jet1_pt = bjet1.Pt();
    b_jet1_eta = bjet1.Eta();
    b_jet2_pt = bjet2.Pt();
    b_jet2_eta = bjet2.Eta();

    if (top1.Pt() < top2.Pt()) { swap(top1, top2); }
    b_top1_pt = top1.Pt();
    b_top1_eta = top1.Eta();
    b_top1_phi = top1.Phi();
    b_top1_rapi = top1.Rapidity();
    b_top1_m = top1.M();
    b_top2_pt = top2.Pt();
    b_top2_eta = top2.Eta();
    b_top2_phi = top2.Phi();
    b_top2_rapi = top2.Rapidity();
    b_top2_m = top2.M();

    auto ttbar = top1+top2;
    b_ttbar_pt = ttbar.Pt();
    b_ttbar_eta = ttbar.Eta();
    b_ttbar_dphi = deltaPhi(top1.Phi(), top2.Phi());
    b_ttbar_m = ttbar.M();
    b_ttbar_rapi = ttbar.Rapidity();


    b_maxweight = maxweight;
    if (maxweight){
      b_step6 = true;
      if (b_step == 5){
        ++b_step;
        cutflow_[9][b_channel]++;
      }
    }
    //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
    // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_met, b_ll_mass, b_maxweight);
    ttree_[sys]->Fill();
  }
}

const reco::Candidate* TtbarDiLeptonAnalyzer::getLast(const reco::Candidate* p) const
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

float TtbarDiLeptonAnalyzer::selectMuons(const cat::MuonCollection& muons, ParticleCollection& selmuons, sys_e sys) const
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
    weight *= mu.scaleFactor("NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1");
    weight *= mu.scaleFactor("NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1");
    selmuons.push_back(mu);
  }
  return weight;
}

float TtbarDiLeptonAnalyzer::selectElecs(const cat::ElectronCollection& elecs, ParticleCollection& selelecs, sys_e sys) const
{
  float weight = 1.;
  for (auto& e : elecs) {
    cat::Electron el(e);
    if (sys == sys_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == sys_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    //if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) continue;
    if ( !el.isTrigMVAValid() or !el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") ) continue;
    if (el.relIso(0.3) > 0.12) continue;

    weight *= el.scaleFactor("mvaEleID-Spring15-25ns-Trig-V1-wp90");
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return weight;
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectJets(const cat::JetCollection& jets, const ParticleCollection& recolep, sys_e sys)
{
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
      if (deltaR(jet.p4(),lep.p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    if (sys == sys_btag_u) b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 1);
    else if (sys == sys_btag_d) b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, -1);
    else b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 0);
    
	seljets.push_back(jet);
  }
  return seljets;
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.605) continue;
    //if (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.89) continue;//forsync
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}
void TtbarDiLeptonAnalyzer::resetBr()
{
  b_nvertex = 0;b_step = -1;b_channel = 0;b_njet = 0;b_nbjet = 0;
  b_step1 = 0;b_step2 = 0;b_step3 = 0;b_step4 = 0;b_step5 = 0;b_step6 = 0;b_tri = 0;b_filtered = 0;
  b_met = -9;
  b_weight = 1; b_puweight = 1; b_puweight_up = 1; b_puweight_dn = 1; b_genweight = 1; b_lepweight = 1;
  b_btagweight = 1;b_btagweight_up = 1;b_btagweight_dn = 1;
  b_pdfWeights.clear();

  b_lep1_pt = -9;b_lep1_eta = -9;b_lep1_phi = -9;
  b_lep2_pt = -9;b_lep2_eta = -9;b_lep2_phi = -9;
  b_ll_pt = -9;b_ll_eta = -9;b_ll_phi = -9;b_ll_m = -9;

  b_partontop1_pt = -9; b_partontop1_eta = -9; b_partontop1_phi = -9; b_partontop1_rapi = -9; b_partontop1_m = -9;
  b_partontop2_pt = -9; b_partontop2_eta = -9; b_partontop2_phi = -9; b_partontop2_rapi = -9; b_partontop2_m = -9;
  b_partonttbar_pt = -9; b_partonttbar_eta = -9; b_partonttbar_dphi = -9; b_partonttbar_m = -9; b_partonttbar_rapi = -9;

  b_partonChannel = -1; b_partonMode1 = -1; b_partonMode2 = -1;
  b_partonlep1_pt = -9; b_partonlep1_eta = -9;
  b_partonlep2_pt = -9; b_partonlep2_eta = -9;
  b_partonjet1_pt = -9; b_partonjet1_eta = -9;
  b_partonjet2_pt = -9; b_partonjet2_eta = -9;
  b_partonInPhase = 0; b_partonInPhaseLep = false; b_partonInPhaseJet = false;

  b_gentop1_pt = -9; b_gentop1_eta = -9; b_gentop1_phi = -9; b_gentop1_rapi = -9; b_gentop1_m = -9;
  b_gentop2_pt = -9; b_gentop2_eta = -9; b_gentop2_phi = -9; b_gentop2_rapi = -9; b_gentop2_m = -9;
  b_genttbar_pt = -9; b_genttbar_eta = -9; b_genttbar_dphi = -9; b_genttbar_m = -9; b_genttbar_rapi = -9;

  b_pseudoTopChannel = -1;
  b_genlep1_pt = -9; b_genlep1_eta = -9;
  b_genlep2_pt = -9; b_genlep2_eta = -9;
  b_genjet1_pt = -9; b_genjet1_eta = -9;
  b_genjet2_pt = -9; b_genjet2_eta = -9;
  b_pseudoInPhase = false;

  b_jet1_pt = -9; b_jet1_eta = -9; b_jet1_CSVInclV2 = -9;
  b_jet2_pt = -9; b_jet2_eta = -9; b_jet2_CSVInclV2 = -9;
  b_top1_pt = -9; b_top1_eta = -9; b_top1_phi = -9; b_top1_rapi = -9; b_top1_m = -9;
  b_top2_pt = -9; b_top2_eta = -9; b_top2_phi = -9; b_top2_rapi = -9; b_top2_m = -9;
  b_ttbar_pt = -9; b_ttbar_eta = -9; b_ttbar_dphi = -9; b_ttbar_m = -9; b_ttbar_rapi = -9;
  b_is3lep = -9;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
