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

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
//#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"
#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstruction.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"
#include "CATTools/CatAnalyzer/interface/analysisUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TTree.h"
#include "TH1D.h"
//#include "TLorentzVector.h"

using namespace std;
using namespace cat;

class TtbarDiLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer();

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

  typedef std::vector<const cat::Lepton*> LeptonPtrs;

  void resetBr();
  float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, sys_e sys) const;
  float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, sys_e sys);
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  const reco::Candidate* getLast(const reco::Candidate* p) const;

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

  BTagWeightEvaluator csvWeight;
  BTagWeightEvaluator bTagWeightL;
  BTagWeightEvaluator bTagWeightM;
  BTagWeightEvaluator bTagWeightT;

  edm::EDGetTokenT<int> recoFiltersToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<vector<float>> pdfweightsToken_, scaleupweightsToken_, scaledownweightsToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightToken_up_, puweightToken_dn_, topPtWeight_;
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
  TH1D * h_nevents;
  int b_run, b_lumi, b_event;
  int b_nvertex, b_step, b_channel, b_njet, b_nbjet;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_step7, b_filtered;
  float b_tri, b_tri_up, b_tri_dn;
  float b_met, b_weight, b_puweight, b_puweight_up, b_puweight_dn, b_genweight,
    b_mueffweight, b_mueffweight_up, b_mueffweight_dn,
    b_eleffweight, b_eleffweight_up, b_eleffweight_dn,
    b_btagweight, b_btagweight_up, b_btagweight_dn;
  float b_topPtWeight;
  //std::vector<float> b_pdfWeights, b_scaleWeights_up, b_scaleWeights_dn, b_csvweights;
  std::vector<float> b_pdfWeights, b_scaleWeights, b_csvweights;
  
  int b_is3lep;
  
  int b_partonChannel, b_partonMode1, b_partonMode2;
  bool b_partonInPhase, b_partonInPhaseLep, b_partonInPhaseJet;
  TLorentzVector b_partonlep1; int b_partonlep1_pid;
  TLorentzVector b_partonlep2; int b_partonlep2_pid;
  TLorentzVector b_partontop1, b_partontop2, b_partonjet1, b_partonjet2, b_partondilep, b_partonttbar;
  float b_partonttbar_dphi;
  
  int b_pseudoChannel;
  bool b_pseudoInPhase;
  TLorentzVector b_pseudolep1; int b_pseudolep1_pid;
  TLorentzVector b_pseudolep2; int b_pseudolep2_pid;
  TLorentzVector b_pseudotop1, b_pseudotop2, b_pseudojet1, b_pseudojet2, b_pseudodilep, b_pseudottbar;
  float b_pseudottbar_dphi;

  TLorentzVector b_lep1; int b_lep1_pid;
  TLorentzVector b_lep2; int b_lep2_pid;
  TLorentzVector b_jet1, b_jet2, b_top1, b_top2, b_dilep, b_ttbar;
  float b_ttbar_dphi;
  float b_jet1_CSVInclV2, b_jet2_CSVInclV2;

  TLorentzVector b_desyjet1, b_desyjet2, b_desytop1, b_desytop2, b_desyttbar;
  float b_desyttbar_dphi;
  float b_desyjet1_CSVInclV2, b_desyjet2_CSVInclV2;

  //std::unique_ptr<TtFullLepKinSolver> solver;
  std::unique_ptr<KinematicSolver> solver_;

  const KinematicReconstruction* kinematicReconstruction;
  typedef math::XYZTLorentzVector LV;
  typedef std::vector<LV> VLV;

  const static int NCutflow = 12;
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
  genWeightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  pdfweightsToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("pdfweights"));	
  scaleupweightsToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaleupweights"));
  scaledownweightsToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaledownweights"));
  topPtWeight_ = consumes<float>(iConfig.getParameter<edm::InputTag>("topPtWeight"));
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
  else if ( algoName == "DEFAULT" ) solver_.reset(new TTDileptonSolver(solverPSet));
  else {
    cerr << "The solver name \"" << solverPSet.getParameter<std::string>("algo") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solver_.reset(new TTDileptonSolver(solverPSet)); // A dummy solver
  }

  csvWeight.initCSVWeight(false, "csvv2");
  bTagWeightL.init(3, "csvv2", BTagEntry::OP_LOOSE , 1);
  bTagWeightM.init(3, "csvv2", BTagEntry::OP_MEDIUM, 1);
  bTagWeightT.init(3, "csvv2", BTagEntry::OP_TIGHT , 1);

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
    tr->Branch("step7", &b_step7, "step7/O");
    tr->Branch("tri", &b_tri, "tri/F");
    tr->Branch("tri_up", &b_tri_up, "tri_up/F");
    tr->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
    tr->Branch("filtered", &b_filtered, "filtered/O");
    tr->Branch("met", &b_met, "met/F");
    tr->Branch("weight", &b_weight, "weight/F");
    tr->Branch("topPtWeight", &b_topPtWeight, "topPtWeight/F");
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
    tr->Branch("btagweight", &b_btagweight, "btagweight/F");
    tr->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
    tr->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");

    if (sys == 0){
      tr->Branch("csvweights","std::vector<float>",&b_csvweights);
      tr->Branch("pdfWeights","std::vector<float>",&b_pdfWeights);
      // tr->Branch("scaleWeights_up","std::vector<float>",&b_scaleWeights_up);
      // tr->Branch("scaleWeights_dn","std::vector<float>",&b_scaleWeights_dn);
      tr->Branch("scaleWeights","std::vector<float>",&b_scaleWeights);
    }
    tr->Branch("is3lep", &b_is3lep, "is3lep/I");

    tr->Branch("partonChannel", &b_partonChannel, "partonChannel/I");
    tr->Branch("partonMode1", &b_partonMode1, "partonMode1/I");    
    tr->Branch("partonMode2", &b_partonMode2, "partonMode2/I");    
    tr->Branch("partonInPhase", &b_partonInPhase, "partonInPhase/O");
    tr->Branch("partonInPhaseLep", &b_partonInPhaseLep, "partonInPhaseLep/O");
    tr->Branch("partonInPhaseJet", &b_partonInPhaseJet, "partonInPhaseJet/O");
    tr->Branch("partonlep1_pid", &b_partonlep1_pid, "partonlep1_pid/I");    
    tr->Branch("partonlep2_pid", &b_partonlep2_pid, "partonlep2_pid/I");
    
    if (sys == 0){
    tr->Branch("partonlep1", "TLorentzVector", &b_partonlep1);
    tr->Branch("partonlep2", "TLorentzVector", &b_partonlep2);
    tr->Branch("partondilep", "TLorentzVector", &b_partondilep);
    tr->Branch("partonjet1", "TLorentzVector", &b_partonjet1);
    tr->Branch("partonjet2", "TLorentzVector", &b_partonjet2);
    tr->Branch("partontop1", "TLorentzVector", &b_partontop1);
    tr->Branch("partontop2", "TLorentzVector", &b_partontop2);
    tr->Branch("partonttbar", "TLorentzVector", &b_partonttbar);
    tr->Branch("partonttbar_dphi", &b_partonttbar_dphi, "partonttbar_dphi/F");
    }    

    tr->Branch("pseudoChannel", &b_pseudoChannel, "pseudoChannel/I");
    tr->Branch("pseudoInPhase", &b_pseudoInPhase, "pseudoInPhase/O");
    tr->Branch("pseudolep1_pid", &b_pseudolep1_pid, "pseudolep1_pid/I");    
    tr->Branch("pseudolep2_pid", &b_pseudolep2_pid, "pseudolep2_pid/I");    

    if (sys == 0){
    tr->Branch("pseudolep1", "TLorentzVector", &b_pseudolep1);
    tr->Branch("pseudolep2", "TLorentzVector", &b_pseudolep2);
    tr->Branch("pseudodilep", "TLorentzVector", &b_pseudodilep);
    tr->Branch("pseudojet1", "TLorentzVector", &b_pseudojet1);
    tr->Branch("pseudojet2", "TLorentzVector", &b_pseudojet2);
    tr->Branch("pseudotop1", "TLorentzVector", &b_pseudotop1);
    tr->Branch("pseudotop2", "TLorentzVector", &b_pseudotop2);
    tr->Branch("pseudottbar", "TLorentzVector", &b_pseudottbar);
    tr->Branch("pseudottbar_dphi", &b_pseudottbar_dphi, "pseudottbar_dphi/F");
    }    

    tr->Branch("lep1", "TLorentzVector", &b_lep1);
    tr->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
    tr->Branch("lep2", "TLorentzVector", &b_lep2);
    tr->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
    tr->Branch("dilep", "TLorentzVector", &b_dilep);
    tr->Branch("jet1", "TLorentzVector", &b_jet1);
    tr->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
    tr->Branch("jet2", "TLorentzVector", &b_jet2);
    tr->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");
    tr->Branch("top1", "TLorentzVector", &b_top1);
    tr->Branch("top2", "TLorentzVector", &b_top2);
    tr->Branch("ttbar", "TLorentzVector", &b_ttbar);
    tr->Branch("ttbar_dphi", &b_ttbar_dphi, "ttbar_dphi/F");

    tr->Branch("desyjet1", "TLorentzVector", &b_desyjet1);
    tr->Branch("desyjet1_CSVInclV2", &b_desyjet1_CSVInclV2, "desyjet1_CSVInclV2/F");
    tr->Branch("desyjet2", "TLorentzVector", &b_desyjet2);
    tr->Branch("desyjet2_CSVInclV2", &b_desyjet2_CSVInclV2, "desyjet2_CSVInclV2/F");
    tr->Branch("desytop1", "TLorentzVector", &b_desytop1);
    tr->Branch("desytop2", "TLorentzVector", &b_desytop2);    
    tr->Branch("desyttbar", "TLorentzVector", &b_desyttbar);
    tr->Branch("desyttbar_dphi", &b_desyttbar_dphi, "desyttbar_dphi/F");
  }

  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});

  kinematicReconstruction = new KinematicReconstruction(1, true);
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
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  const bool runOnMC = !iEvent.isRealData();
  cutflow_[0][0]++;

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();
    bool keepTtbarSignal = false;

    edm::Handle<int> partonTop_channel;
    if ( iEvent.getByToken(partonTop_channel_, partonTop_channel)){

      edm::Handle<float> topPtWeightHandle;
      iEvent.getByToken(topPtWeight_, topPtWeightHandle);
      b_topPtWeight = *topPtWeightHandle;

      if (sys == sys_nom){
        edm::Handle<vector<float>> pdfweightsHandle;
        iEvent.getByToken(pdfweightsToken_, pdfweightsHandle);
        for (const float & aPdfWeight : *pdfweightsHandle){
          b_pdfWeights.push_back(aPdfWeight);
        }
        edm::Handle<vector<float>> scaleupweightsHandle;
        iEvent.getByToken(scaleupweightsToken_, scaleupweightsHandle);
        for (const float & aScaleWeight : *scaleupweightsHandle){
          b_scaleWeights.push_back(aScaleWeight);
        }
        edm::Handle<vector<float>> scaledownweightsHandle;
        iEvent.getByToken(scaledownweightsToken_, scaledownweightsHandle);
        for (const float & aScaleWeight : *scaledownweightsHandle){
          b_scaleWeights.push_back(aScaleWeight);
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
      if (b_partonChannel == CH_FULLLEPTON) keepTtbarSignal = true;

      if ( !(partonTop_genParticles->empty()) ){

        // Get Top quark pairs
	auto parton1 = &partonTop_genParticles->at(0);
	auto parton2 = &partonTop_genParticles->at(1);
	if (parton1->charge() < 0) swap(parton1, parton2);
	
        b_partontop1 = ToTLorentzVector(*parton1);
        b_partontop2 = ToTLorentzVector(*parton2);
	b_partonttbar = b_partontop1 + b_partontop2;
	b_partonttbar_dphi = b_partontop1.DeltaPhi(b_partontop2);

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
	    b_partonlep1 = ToTLorentzVector(*partonW11);
	    b_partonlep2 = ToTLorentzVector(*partonW21);
	    b_partonlep1_pid = partonW11->pdgId();
	    b_partonlep2_pid = partonW21->pdgId();
	    b_partondilep = b_partonlep1 + b_partonlep2;
	    b_partonjet1 = ToTLorentzVector(*partonB1);
	    b_partonjet2 = ToTLorentzVector(*partonB2);
          }
        }
        if (b_partonInPhaseJet && b_partonInPhaseLep) b_partonInPhase = true;
      }

      // Start to build pseudo top
      b_pseudoChannel = CH_NOLL;
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
          case 22: b_pseudoChannel = CH_ELEL; break;
          case 26: b_pseudoChannel = CH_MUMU; break;
          case 24: b_pseudoChannel = CH_MUEL; break;
          default: b_pseudoChannel = CH_NOLL;
        }
	if (b_pseudoChannel > 0) keepTtbarSignal = true;
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
	
        // if ( true ) {
        //   const auto t1Alt = w1 + bjet2;
        //   const auto t2Alt = w2 + bjet1;

        //   const double tMass = 172.5;
        //   const double dm = std::abs(gentop1.mass()-tMass)+std::abs(gentop2.mass()-tMass);
        //   const double dmAlt = std::abs(t1Alt.mass()-tMass)+std::abs(t2Alt.mass()-tMass);
        //   if ( dm > dmAlt ) { gentop1 = t1Alt; gentop2 = t2Alt; std::swap(bjet1, bjet2); }
        // }
	//        if (gentop1.Pt() < gentop2.Pt()) { swap(gentop1, gentop2); }
	if (pseudoTopLeptonHandle->at(leptonIdxs[0]).charge() < 0) swap(gentop1, gentop2);
        b_pseudotop1 = ToTLorentzVector(gentop1);
        b_pseudotop2 = ToTLorentzVector(gentop2);
	b_pseudottbar = b_pseudotop1 + b_pseudotop2;
	b_pseudottbar_dphi = b_pseudotop1.DeltaPhi(b_pseudotop2);
	
        b_pseudolep1 = ToTLorentzVector(pseudoTopLeptonHandle->at(leptonIdxs[0]));
        b_pseudolep2 = ToTLorentzVector(pseudoTopLeptonHandle->at(leptonIdxs[1]));
	b_pseudodilep = b_pseudolep1 + b_pseudolep2;
	//b_pseudolep1_pid = lepton1.pdgId();
	//b_pseudolep2_pid = lepton2.pdgId();
	
        if (bjet1.Pt() < bjet2.Pt()) { swap(bjet1, bjet2); }
        b_pseudojet1 = ToTLorentzVector(bjet1);
        b_pseudojet2 = ToTLorentzVector(bjet2);

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
      iEvent.getByToken(genWeightToken_, genweightHandle);
      b_genweight = (*genweightHandle);
      b_weight = b_genweight*b_puweight;
    }

    if (sys == sys_nom) h_nevents->Fill(0.5,b_puweight*b_genweight);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()){ // skip the event if no PV found
      if (keepTtbarSignal) ttree_[sys]->Fill();
      continue;
    }
    if (sys == sys_nom) cutflow_[1][b_channel]++;

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
      if (keepTtbarSignal) ttree_[sys]->Fill();
      continue;
    }
    if (sys == sys_nom) cutflow_[3][b_channel]++;

    std::vector<const cat::Lepton*> recolep;
    for ( const auto& x : selMuons ) recolep.push_back(&x);
    for ( const auto& x : selElecs ) recolep.push_back(&x);

    sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
    recolep.erase(recolep.begin()+2,recolep.end());
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
    b_tri = b_tri_up = b_tri_dn = 0;
    edm::Handle<int> trigHandle;
    if      ( b_channel == CH_ELEL ) iEvent.getByToken(trigTokenELEL_, trigHandle);
    else if ( b_channel == CH_MUMU ) iEvent.getByToken(trigTokenMUMU_, trigHandle);
    else if ( b_channel == CH_MUEL ) iEvent.getByToken(trigTokenMUEL_, trigHandle);
    if ( *trigHandle != 0 ) {
       b_tri = computeTrigSF(recolep1, recolep2);
       b_tri_up = computeTrigSF(recolep1, recolep2,  1);
       b_tri_dn = computeTrigSF(recolep1, recolep2, -1);
    }

    b_lep1 = recolep1.tlv(); b_lep1_pid = recolep1.pdgId();
    b_lep2 = recolep2.tlv(); b_lep2_pid = recolep2.pdgId();
    b_dilep = b_lep1+b_lep2;
    const auto tlv_ll = recolep1.p4()+recolep2.p4();

    if (tlv_ll.M() < 20. || recolep1.charge() * recolep2.charge() > 0){
      if (keepTtbarSignal) ttree_[sys]->Fill();
      continue;
    }
    b_step1 = true;
    b_step = 1;
    if (sys == sys_nom) cutflow_[4][b_channel]++;

    if ( (b_channel == CH_MUEL) || ((tlv_ll.M() < 76) || (tlv_ll.M() > 106)) ){
      b_step2 = true;
      b_step = 2;
      if (sys == sys_nom) cutflow_[5][b_channel]++;
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
	
	b_step7 = true;	
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
    ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
    //int kin=0;
    //const cat::Jet* kinj1, * kinj2;

    const auto recolepLV1= recolep1.p4();
    const auto recolepLV2= recolep2.p4();
    std::vector<cat::KinematicSolution> sol2Bs, sol1Bs, sol0Bs;
    cat::KinematicSolution bestSol;

    for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
      const auto recojet1= jet1->p4();
      const bool isBjet1 = jet1->bDiscriminator(BTAG_CSVv2) >= WP_BTAG_CSVv2L;
      for (auto jet2 = next(jet1); jet2 != end; ++jet2){
        const auto recojet2= jet2->p4();
        const bool isBjet2 = jet2->bDiscriminator(BTAG_CSVv2) >= WP_BTAG_CSVv2L;

        solver_->solve(met, recolepLV1, recolepLV2, recojet2, recojet1);
        const cat::KinematicSolution sol1 = solver_->solution();

        solver_->solve(met, recolepLV1, recolepLV2, recojet1, recojet2);
        const cat::KinematicSolution sol2 = solver_->solution();

        //if      ( sol1.quality() >= sol2.quality() and sol1.quality() > bestSol.quality() ) bestSol = sol1;
        //else if ( sol2.quality() >= sol1.quality() and sol2.quality() > bestSol.quality() ) bestSol = sol2;

        if ( isBjet1 and isBjet2 ) {
          sol2Bs.push_back(sol1);
          sol2Bs.push_back(sol2);
        }
        else if ( isBjet1 or isBjet2 ) {
          sol1Bs.push_back(sol1);
          sol1Bs.push_back(sol2);
        }
        else {
          sol0Bs.push_back(sol1);
          sol0Bs.push_back(sol2);
        }
      }
    }
    auto greaterByQuality = [](const cat::KinematicSolution& a, const cat::KinematicSolution& b) { return a.quality() > b.quality(); };
    std::sort(sol2Bs.begin(), sol2Bs.end(), greaterByQuality);
    std::sort(sol1Bs.begin(), sol1Bs.end(), greaterByQuality);
    std::sort(sol0Bs.begin(), sol0Bs.end(), greaterByQuality);

    // Prefer to select b jets
    if      ( !sol2Bs.empty() ) bestSol = sol2Bs.front();
    else if ( !sol1Bs.empty() ) bestSol = sol1Bs.front();
    //else if ( !sol0Bs.empty() ) bestSol = sol0Bs.front();

    //saving results
    const double maxweight = bestSol.quality();
    if ( maxweight > 0 ) {
      math::XYZTLorentzVector top1, top2, nu1, nu2, bjet1, bjet2;
      bjet1 = bestSol.j1();
      bjet2 = bestSol.j2();
      top1 = bestSol.t1();
      top2 = bestSol.t2();
      if (recolep1.charge() < 0) swap(top1, top2);

      if (bjet1.Pt() < bjet2.Pt()) { swap(bjet1, bjet2); }


      b_jet1 = ToTLorentzVector(bjet1);
      b_jet2 = ToTLorentzVector(bjet2);
      b_top1 = ToTLorentzVector(top1);
      b_top2 = ToTLorentzVector(top2);
      b_ttbar = b_top1 + b_top2;
      b_ttbar_dphi = b_top1.DeltaPhi(b_top2);

      b_step6 = true;
      if (b_step == 5){
        ++b_step;
        if (sys == sys_nom) cutflow_[9][b_channel]++;
      }
      //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
      // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_met, b_ll_mass, b_maxweight);
    }
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

float TtbarDiLeptonAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, sys_e sys) const
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

float TtbarDiLeptonAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
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

cat::JetCollection TtbarDiLeptonAnalyzer::selectJets(const cat::JetCollection& jets, const TtbarDiLeptonAnalyzer::LeptonPtrs& recolep, sys_e sys)
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
    //if (sys == sys_btag_u) b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 1);
    //else if (sys == sys_btag_d) b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, -1);
    //else b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 0);
    for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] *= csvWeight.getSF(jet, iu);
    seljets.push_back(jet);
  }
  for (unsigned int iu=0; iu<19; iu++) b_csvweights.push_back(Jet_SF_CSV[iu]);

  b_btagweight = Jet_SF_CSV[0];
  // if      ( sys == sys_btag_u ) b_btagweight = bTagWeightL.eventWeight(seljets, 1);
  // else if ( sys == sys_btag_d ) b_btagweight = bTagWeightL.eventWeight(seljets, 2);
  // else                          b_btagweight = bTagWeightL.eventWeight(seljets, 0);

  return seljets;
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectBJets(const JetCollection& jets) const
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
void TtbarDiLeptonAnalyzer::resetBr()
{
  b_channel = 0;
  b_nvertex = b_njet = b_nbjet = 0;
  b_step = -1;
  b_step1 = b_step2 = b_step3 = b_step4 = b_step5 = b_step6 = b_step7 = 0;
  b_tri = b_tri_up = b_tri_dn = 0;
  b_filtered = 0;
  b_met = -9;
  b_weight = b_puweight = b_puweight_up = b_puweight_dn = b_genweight = 1;
  b_mueffweight = b_mueffweight_up = b_mueffweight_dn = 1;
  b_eleffweight = b_eleffweight_up = b_eleffweight_dn = 1;
  b_btagweight = b_btagweight_up = b_btagweight_dn = 1;
  b_topPtWeight = 1.;
  b_csvweights.clear();
  b_pdfWeights.clear();
  b_scaleWeights.clear();
  //b_scaleWeights_up.clear(); b_scaleWeights_dn.clear();

  b_is3lep = -9;
  
  b_partonChannel = 0; b_partonMode1 = 0; b_partonMode2 = 0;
  b_partonInPhase = 0; b_partonInPhaseLep = 0; b_partonInPhaseJet = 0;
  b_partonlep1 = TLorentzVector(); b_partonlep1_pid = 0;
  b_partonlep2 = TLorentzVector(); b_partonlep2_pid = 0;
  b_partontop1 = TLorentzVector(); b_partontop2 = TLorentzVector();
  b_partonjet1 = TLorentzVector(); b_partonjet2 = TLorentzVector();
  b_partondilep = TLorentzVector(); b_partonttbar = TLorentzVector();
  b_partonttbar_dphi = 0;
  
  b_pseudoChannel = 0;
  b_pseudoInPhase = 0;
  b_pseudolep1 = TLorentzVector(); b_pseudolep1_pid = 0;
  b_pseudolep2 = TLorentzVector(); b_pseudolep2_pid = 0;
  b_pseudotop1 = TLorentzVector(); b_pseudotop2 = TLorentzVector();
  b_pseudojet1 = TLorentzVector(); b_pseudojet2 = TLorentzVector();
  b_pseudodilep = TLorentzVector(); b_pseudottbar = TLorentzVector();
  b_pseudottbar_dphi = 0;

  b_lep1 = TLorentzVector(); b_lep1_pid = 0;
  b_lep2 = TLorentzVector(); b_lep2_pid = 0;
  b_jet1 = TLorentzVector(); b_jet2 = TLorentzVector();
  b_top1 = TLorentzVector(); b_top2 = TLorentzVector();
  b_jet1_CSVInclV2 = 0; b_jet2_CSVInclV2 = 0;
  b_dilep = TLorentzVector(); b_ttbar = TLorentzVector();
  b_ttbar_dphi = 0;

  b_desyjet1 = TLorentzVector(); b_desyjet2 = TLorentzVector();
  b_desytop1 = TLorentzVector(); b_desytop2 = TLorentzVector();
  b_desyjet1_CSVInclV2 = 0; b_desyjet2_CSVInclV2 = 0;
  b_desyttbar = TLorentzVector();
  b_desyttbar_dphi = 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
