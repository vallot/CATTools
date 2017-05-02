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
//#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
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

//########### Dstar begin
#include "CATTools/DataFormats/interface/SecVertex.h"
#include<TClonesArray.h>
#include"dileptonCommon.h"


using namespace std;
using namespace cat;
using namespace dileptonCommonGlobal;

//########### Dstar end

class TtbarDiLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  void setBranch(TTree* tree, systematic sys);
  void resetBranch();
  bool genInformation(const edm::Event& iEvent);
  bool eventSelection(const edm::Event& iEvent, systematic sys);
  
  cat::MuonCollection selectMuons(const cat::MuonCollection& muons, systematic sys);
  cat::ElectronCollection selectElecs(const cat::ElectronCollection& elecs, systematic sys);
  typedef std::vector<const cat::Lepton*> LeptonPtrs;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, systematic sys);
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  const reco::Candidate* getLast(const reco::Candidate* p) const;
  int isFromtop( const reco::GenParticle& p);
  shared_ptr<TLorentzVector> mcMatching( vector<TLorentzVector>& , TLorentzVector& ) ;

  ScaleFactorEvaluator muonSF_, elecSF_;

  BTagWeightEvaluator csvWeight, bTagWeightL, bTagWeightM, bTagWeightT;

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

  bool runOnMC;
  
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

  // ####### Dstar  begin  #########

  std::vector<bool> b_d0_true, b_d0_fit;
  //std::vector<float> b_d0_pt, b_d0_eta, b_d0_phi, b_d0_m;
  std::vector<float> b_d0_LXY, b_d0_L3D, b_d0_dRTrue, b_d0_relPtTrue, b_d0_dca;
  std::vector<float> b_d0_dau1_q;
  std::vector<float> b_d0_dau2_q;
  std::vector<float> b_d0_vProb;

  std::vector<float> b_d0_lepSV_lowM;
  std::vector<float> b_d0_lepSV_dRM;
  std::vector<float> b_d0_lepSV_correctM; // for test

  std::vector<bool> b_dstar_true, b_dstar_fit;
  //std::vector<float> b_dstar_pt, b_dstar_eta, b_dstar_phi, b_dstar_m;
  std::vector<float> b_dstar_q, b_dstar_LXY, b_dstar_L3D, b_dstar_dRTrue, b_dstar_relPtTrue, b_dstar_dca, b_dstar_dca2, b_dstar_dca3;
  std::vector<float> b_dstar_dau1_q;
  std::vector<float> b_dstar_dau2_q;
  std::vector<float> b_dstar_dau3_q;
  std::vector<float> b_dstar_vProb;
  std::vector<float> b_dstar_diffMass;

  std::vector<float> b_dstar_lepSV_lowM;
  std::vector<float> b_dstar_lepSV_dRM;
  std::vector<float> b_dstar_opCharge_M;
  std::vector<float> b_dstar_lepSV_correctM;

  std::vector<bool> b_Jpsi_true, b_Jpsi_fit;
  //std::vector<float> b_Jpsi_pt, b_Jpsi_eta, b_Jpsi_phi, b_Jpsi_m;
  std::vector<float> b_Jpsi_LXY, b_Jpsi_L3D, b_Jpsi_dRTrue, b_Jpsi_relPtTrue, b_Jpsi_dca;
  std::vector<float> b_Jpsi_dau1_q;
  std::vector<float> b_Jpsi_dau2_q;
  std::vector<int>   b_Jpsi_dau_pid;
  std::vector<float> b_Jpsi_vProb;

  std::vector<float> b_Jpsi_lepSV_lowM;
  std::vector<float> b_Jpsi_lepSV_dRM;
  std::vector<float> b_Jpsi_lepSV_correctM; // for test

  TClonesArray *b_d0,    *b_d0_dau1,    *b_d0_dau2;
  TClonesArray *b_dstar, *b_dstar_dau1, *b_dstar_dau2, *b_dstar_dau3;
  TClonesArray *b_Jpsi,    *b_Jpsi_dau1,    *b_Jpsi_dau2;

  edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
  edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;
  edm::EDGetTokenT<cat::SecVertexCollection>      JpsiToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >            mcSrc_;

  double matchingDeltaR_;
  // ############### Dstar end #################
};


TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
// ############## Dstar begin #####################
  //parameterInit(iConfig);

  d0Token_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("d0s"));
  dstarToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("dstars"));
  JpsiToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("Jpsis"));
  mcSrc_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("mcLabel"));
  matchingDeltaR_  = iConfig.getParameter<double>("matchingDeltaR");


//  for (int sys = 0; sys < nsys_e; ++sys){
//    auto tr = ttree_[sys];
//    setBranch(tr, sys);
//  }


// ############## Dstar end #####################
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
              muonSFSet.getParameter<vdouble>("eta_bins"),
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
  h_nevents = fs->make<TH1D>("nevents","nevents",1,0,1);       
  for (int sys = 0; sys < syst_total; ++sys){
    ttree_.push_back(fs->make<TTree>(systematicName[systematic(sys)].c_str(), systematicName[systematic(sys)].c_str()));
    auto tr = ttree_.back();
    setBranch(tr, systematic(sys));
  }
  
  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});

  kinematicReconstruction = new KinematicReconstruction(1, true);

}

TtbarDiLeptonAnalyzer::~TtbarDiLeptonAnalyzer()
{
  cout <<"      cut flow   emu    ee    mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step "<< i << "    "<< cutflow_[i][0] <<  "   "<< cutflow_[i][1] << "   " << cutflow_[i][2] << "   " << cutflow_[i][3]<< endl;
  }
  cout << "          step1   step2      step3      step4     step5" << endl;
  cout << "MuEl" <<"        "<< cutflow_[4][1] <<"    "<< cutflow_[5][1]<<"       "<<cutflow_[6][1]<<"        "<<cutflow_[7][1]<<"        "<<cutflow_[8][1]<<endl;
  cout << "ElEl" <<"        "<< cutflow_[4][2] <<"    "<< cutflow_[5][2]<<"       "<<cutflow_[6][2]<<"        "<<cutflow_[7][2]<<"        "<<cutflow_[8][2]<<endl;
  cout << "MuMu" <<"        "<< cutflow_[4][3] <<"    "<< cutflow_[5][3]<<"       "<<cutflow_[6][3]<<"        "<<cutflow_[7][3]<<"        "<<cutflow_[8][3]<<endl;
}


// #####################  Dstar begin  ########################

int TtbarDiLeptonAnalyzer::isFromtop( const reco::GenParticle& p){
  int nIDMother;

  int nTopMother = 0;

  //  string pt = Form("%f", p.pt());
  //  string pdgid = Form("%i",p.pdgId());
  const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(p.mother());
  while( mother != nullptr ) {
    //    string id = Form("%i", mother->pdgId());
    //    string mopt = Form("%f", mother->pt());
    nIDMother = mother->pdgId();

    if( abs(nIDMother) == 6 ) {
      nTopMother = nIDMother;
      break;
    }

    mother = dynamic_cast<const reco::GenParticle*>(mother->mother());
  }

  return nTopMother;
}

shared_ptr<TLorentzVector> TtbarDiLeptonAnalyzer::mcMatching( vector<TLorentzVector>& aGens, TLorentzVector& aReco) {
  float minDR= 999.;
  //float minRelPt = 1.0;
  shared_ptr<TLorentzVector> matchedGen;
  for( auto& aGen : aGens ) {
    float deltaR = aGen.DeltaR(aReco);
    if ( deltaR < minDR ) { matchedGen=make_shared<TLorentzVector>(aGen); minDR = deltaR; }
  }
  if ( minDR < matchingDeltaR_ ) {
    return matchedGen;
  }
  return nullptr;
}



// ###################### Dstar end #########################



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
  runOnMC = !iEvent.isRealData();

  bool keepGenSignal = false;
  bool keepEvent = false;
  
  for (int sys = 0; sys < syst_total; ++sys){
    if (sys != syst_nom && !runOnMC) break;
    resetBranch();
    
    if (runOnMC && sys == syst_nom)
      keepGenSignal = genInformation(iEvent);
    
    keepEvent = eventSelection(iEvent, systematic(sys));
    
    if (keepGenSignal || keepEvent)
      ttree_[sys]->Fill();
  
  }
}

bool TtbarDiLeptonAnalyzer::eventSelection(const edm::Event& iEvent, systematic sys)
{
  if (sys == syst_nom) cutflow_[0][0]++;
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
    iEvent.getByToken(genWeightToken_, genweightHandle);
    b_genweight = (*genweightHandle);
    b_weight = b_genweight*b_puweight;
    
    h_nevents->Fill(0.5,b_weight);
  }  
// #####################  Dstar begin  ########################
  edm::Handle<cat::SecVertexCollection> d0s;       iEvent.getByToken(d0Token_,d0s);
  edm::Handle<cat::SecVertexCollection> dstars;    iEvent.getByToken(dstarToken_,dstars);
  edm::Handle<cat::SecVertexCollection> Jpsis;     iEvent.getByToken(JpsiToken_,Jpsis);

  edm::Handle<edm::View<reco::GenParticle> > mcHandle;

  if ( runOnMC ) {
    iEvent.getByToken(mcSrc_, mcHandle);
  }

  vector<TLorentzVector> gen_d0s;
  vector<TLorentzVector> gen_dstars;
  vector<TLorentzVector> gen_Jpsis;

  int nIDMother;

  TLorentzVector vecSumMom;

  if ( runOnMC ) {
    for( const auto& aGenParticle : *mcHandle) {
      //If genParticle is D0,
      if ( std::abs(aGenParticle.pdgId()) == 421 ) {
        gen_d0s.push_back( ToTLorentzVector(aGenParticle));

        nIDMother = isFromtop(aGenParticle);

        if ( /*0 == 1 &&*/ abs(nIDMother) != 6 ) {
            continue;
        }

        //printf("Sign : %i\n", nIDMother * aGenParticle.pdgId() / abs(aGenParticle.pdgId()));
        //printf("%s\n", ( aGenParticle.charge() * aGenParticle.pdgId() >= 0 ? "S" : "O" ));
        vecSumMom = ToTLorentzVector(aGenParticle) +
          ( nIDMother * b_lep1_pid < 0 ? b_lep1 : b_lep2 );
        b_d0_lepSV_correctM.push_back(vecSumMom.M());
      } else if ( std::abs(aGenParticle.pdgId()) ==  413 ) {
        gen_dstars.push_back( ToTLorentzVector(aGenParticle));

        nIDMother = isFromtop(aGenParticle);

        if ( /*0 == 1 &&*/ abs(nIDMother) != 6 ) {
            continue;
        }

        //printf("Sign : %i\n", nIDMother * aGenParticle.pdgId() / abs(aGenParticle.pdgId()));
        //printf("%s\n", ( aGenParticle.charge() * aGenParticle.pdgId() >= 0 ? "S" : "O" ));

        vecSumMom = ToTLorentzVector(aGenParticle) +
          ( nIDMother * b_lep1_pid < 0 ? b_lep1 : b_lep2 );
        b_dstar_lepSV_correctM.push_back(vecSumMom.M());
      } else if ( std::abs(aGenParticle.pdgId()) ==  443 ) {
        gen_Jpsis.push_back( ToTLorentzVector(aGenParticle));

        nIDMother = isFromtop(aGenParticle);

        if ( /*0 == 1 &&*/ abs(nIDMother) != 6 ) {
            continue;
        }

        //printf("Sign : %i\n", nIDMother * aGenParticle.pdgId() / abs(aGenParticle.pdgId()));
        //printf("%s\n", ( aGenParticle.charge() * aGenParticle.pdgId() >= 0 ? "S" : "O" ));

        vecSumMom = ToTLorentzVector(aGenParticle) +
          ( nIDMother * b_lep1_pid < 0 ? b_lep1 : b_lep2 );
        b_Jpsi_lepSV_correctM.push_back(vecSumMom.M());
      }
    }
  }
  int d0_count=-1;
  int dstar_count=-1;
  int Jpsi_count=-1;

  TClonesArray& br_d0 = *b_d0;
  TClonesArray& br_d0_dau1 = *b_d0_dau1;
  TClonesArray& br_d0_dau2 = *b_d0_dau2;

  TClonesArray& br_dstar = *b_dstar;
  TClonesArray& br_dstar_dau1 = *b_dstar_dau1;
  TClonesArray& br_dstar_dau2 = *b_dstar_dau2;
  TClonesArray& br_dstar_dau3 = *b_dstar_dau3;

  TClonesArray& br_Jpsi = *b_Jpsi;
  TClonesArray& br_Jpsi_dau1 = *b_Jpsi_dau1;
  TClonesArray& br_Jpsi_dau2 = *b_Jpsi_dau2;

  TLorentzVector vecDMMom, vecDau12;
  float fQDau, fQDM;

  TLorentzVector vecSumDMLep1, vecSumDMLep2;
  float fMDMLep1, fMDMLep2;

  float fDeltaEta, fDeltaPhi;
  float fSqrtdRMLep1, fSqrtdRMLep2;

  for( auto& x : *d0s) {
    d0_count++;

    auto d0_tlv = ToTLorentzVector(x);
    new( br_d0[d0_count]) TLorentzVector(d0_tlv);
    new( br_d0_dau1[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0))));
    new( br_d0_dau2[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));

    b_d0_dca.push_back( x.dca());
    double d0_vProb = x.vProb();
    b_d0_vProb.push_back(d0_vProb);

    if ( abs( d0_vProb ) > 1e-5 ) {
      b_d0_fit.push_back(true);
      b_d0_L3D.push_back( x.l3D());
      b_d0_LXY.push_back( x.lxy());
    } else {
      b_d0_fit.push_back(false);
      b_d0_L3D.push_back( -9 );
      b_d0_LXY.push_back( -9 );
    }

    if ( runOnMC ) {
        shared_ptr<TLorentzVector> genMatched = mcMatching( gen_d0s, d0_tlv );
        if ( genMatched != nullptr) {
          b_d0_true.push_back( true );
          b_d0_dRTrue.push_back( genMatched->DeltaR( d0_tlv ));
          b_d0_relPtTrue.push_back( (genMatched->Pt()- d0_tlv.Pt())/genMatched->Pt());
        }
        else {
          b_d0_true.push_back( false );
          b_d0_dRTrue.push_back( -9);
          b_d0_relPtTrue.push_back(-9);
        }
    }

    fQDM = 0;

    fQDau = x.daughter(0)->charge();
    b_dstar_dau1_q.push_back(fQDau);
    fQDM += fQDau;
    fQDau = x.daughter(0)->charge();
    b_dstar_dau1_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(1)->charge();
    b_dstar_dau2_q.push_back(fQDau);
    fQDM += fQDau;

    vecDMMom = ToTLorentzVector(*(x.daughter(0))) + ToTLorentzVector(*(x.daughter(1)));

    vecSumDMLep1 = b_lep1 + vecDMMom;
    vecSumDMLep2 = b_lep2 + vecDMMom;

    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();

    fDeltaEta = b_lep1.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep1.Phi() - vecDMMom.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    fDeltaEta = b_lep2.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep2.Phi() - vecDMMom.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    b_d0_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_d0_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
  }
  for( auto& x : *dstars) {
    dstar_count++;

    auto dstar_tlv = ToTLorentzVector(x);
    new( br_dstar[dstar_count]) TLorentzVector( dstar_tlv );
    new( br_dstar_dau1[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0))));
    new( br_dstar_dau2[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));
    new( br_dstar_dau3[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(2))));

    b_dstar_dca.push_back( x.dca());
    b_dstar_dca2.push_back( x.dca(1));
    b_dstar_dca3.push_back( x.dca(2));

    double dstar_vProb = x.vProb();
    b_dstar_vProb.push_back(dstar_vProb);

    if ( abs( dstar_vProb) > 1e-5) {
      b_dstar_fit.push_back(true);
      b_dstar_L3D.push_back( x.l3D());
      b_dstar_LXY.push_back( x.lxy());
    } else {
      b_dstar_fit.push_back(false);
      b_dstar_L3D.push_back( -9 );
      b_dstar_LXY.push_back( -9 );
    }
    if ( runOnMC ) {
        shared_ptr<TLorentzVector> genMatched = mcMatching( gen_dstars, dstar_tlv );
        if ( genMatched != nullptr) {
          b_dstar_true.push_back( true );
          b_dstar_dRTrue.push_back( genMatched->DeltaR( dstar_tlv));
          b_dstar_relPtTrue.push_back( (genMatched->Pt()- dstar_tlv.Pt())/genMatched->Pt());
        }
        else {
          b_dstar_true.push_back( false );
          b_dstar_dRTrue.push_back( -9);
          b_dstar_relPtTrue.push_back(-9);
        }
    }

    fQDM = 0;

    fQDau = x.daughter(0)->charge();
    b_dstar_dau1_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(1)->charge();
    b_dstar_dau2_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(2)->charge();
    b_dstar_dau3_q.push_back(fQDau);
    fQDM += fQDau;

    vecDMMom = ToTLorentzVector(*(x.daughter(0))) +
        ToTLorentzVector(*(x.daughter(1))) +
        ToTLorentzVector(*(x.daughter(2)));

    vecDau12 = ToTLorentzVector(*(x.daughter(0))) +
        ToTLorentzVector(*(x.daughter(1)));
    b_dstar_diffMass.push_back(vecDMMom.M() - vecDau12.M());

    vecSumDMLep1 = b_lep1 + vecDMMom;
    vecSumDMLep2 = b_lep2 + vecDMMom;

    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();

    fDeltaEta = b_lep1.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep1.Phi() - vecDMMom.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    fDeltaEta = b_lep2.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep2.Phi() - vecDMMom.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    b_dstar_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_dstar_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_dstar_opCharge_M.push_back(( fQDM * b_lep1_pid <= 0.0 ? fMDMLep1 : fMDMLep2 ));
  }
  for( auto& x : *Jpsis) {
    Jpsi_count++;

    auto Jpsi_tlv = ToTLorentzVector(x);
    new( br_Jpsi[Jpsi_count]) TLorentzVector( Jpsi_tlv );
    new( br_Jpsi_dau1[Jpsi_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0))));
    new( br_Jpsi_dau2[Jpsi_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));

    b_Jpsi_dca.push_back( x.dca());

    double Jpsi_vProb = x.vProb();
    b_Jpsi_vProb.push_back(Jpsi_vProb);
    if ( abs( Jpsi_vProb) > 1e-5) {
      b_Jpsi_fit.push_back(true);
      b_Jpsi_L3D.push_back( x.l3D());
      b_Jpsi_LXY.push_back( x.lxy());
    } else {
      b_Jpsi_fit.push_back(false);
      b_Jpsi_L3D.push_back( -9 );
      b_Jpsi_LXY.push_back( -9 );
    }

    if ( runOnMC ) {
        shared_ptr<TLorentzVector> genMatched = mcMatching( gen_Jpsis, Jpsi_tlv );
        if ( genMatched != nullptr) {
          b_Jpsi_true.push_back( true );
          b_Jpsi_dRTrue.push_back( genMatched->DeltaR( Jpsi_tlv));
          b_Jpsi_relPtTrue.push_back( (genMatched->Pt()- Jpsi_tlv.Pt())/genMatched->Pt());
        }
        else {
          b_Jpsi_true.push_back( false );
          b_Jpsi_dRTrue.push_back( -9);
          b_Jpsi_relPtTrue.push_back(-9);
        }
    }

    fQDM = 0;

    fQDau = x.daughter(0)->charge();
    b_Jpsi_dau1_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(1)->charge();
    b_Jpsi_dau2_q.push_back(fQDau);
    fQDM += fQDau;

    b_Jpsi_dau_pid.push_back(x.daughter(0)->pdgId());

    vecDMMom = ToTLorentzVector(*(x.daughter(0))) + ToTLorentzVector(*(x.daughter(1)));
    vecSumDMLep1 = b_lep1 + vecDMMom;
    vecSumDMLep2 = b_lep2 + vecDMMom;

    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();

    fDeltaEta = b_lep1.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep1.Phi() - vecDMMom.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    fDeltaEta = b_lep2.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep2.Phi() - vecDMMom.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;

    b_Jpsi_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_Jpsi_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
  }
  

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) // skip the event if no PV found
    return false;
  
  if (sys == syst_nom) cutflow_[1][b_channel]++;

  // const reco::Vertex &PV = vertices->front();
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  b_nvertex = *nGoodVertexHandle;

  edm::Handle<int> lumiSelectionHandle;
  iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
  if (!runOnMC){
    if (*lumiSelectionHandle == 0) return false;
  }

  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;

  if (sys == syst_nom) cutflow_[2][b_channel]++;

  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);

  // Find leptons and sort by pT
  cat::MuonCollection selMuons = selectMuons(*muons, sys);
  cat::ElectronCollection selElecs = selectElecs(*electrons, sys);
    
  if ( selMuons.size()+selElecs.size() < 2 ) return false;
  
  if (sys == syst_nom) cutflow_[3][b_channel]++;

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

  b_mueffweight    = muonSF_.getScaleFactor(recolep1, 13, 0)*muonSF_.getScaleFactor(recolep2, 13,  0);
  b_mueffweight_up = muonSF_.getScaleFactor(recolep1, 13, +1)*muonSF_.getScaleFactor(recolep2, 13, +1);
  b_mueffweight_dn = muonSF_.getScaleFactor(recolep1, 13, -1)*muonSF_.getScaleFactor(recolep2, 13, -1);

  b_eleffweight    = elecSF_.getScaleFactor(recolep1, 11, 0)*elecSF_.getScaleFactor(recolep2, 11,  0);
  b_eleffweight_up = elecSF_.getScaleFactor(recolep1, 11, +1)*elecSF_.getScaleFactor(recolep2, 11, +1);
  b_eleffweight_dn = elecSF_.getScaleFactor(recolep1, 11, -1)*elecSF_.getScaleFactor(recolep2, 11, -1);

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

  if (tlv_ll.M() < 20. || recolep1.charge() * recolep2.charge() > 0) return false;
  
  b_step1 = true;
  b_step = 1;
  if (sys == syst_nom) cutflow_[4][b_channel]++;

  if ( (b_channel == CH_MUEL) || ((tlv_ll.M() < 76) || (tlv_ll.M() > 106)) ){
    b_step2 = true;
    b_step = 2;
    if (sys == syst_nom) cutflow_[5][b_channel]++;
  }

  JetCollection&& selectedJets = selectJets(*jets, recolep, sys);
  JetCollection&& selectedBJets = selectBJets(selectedJets);

  const auto met = mets->front().p4();
  b_met = met.pt();
  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();

  if ((b_channel == CH_MUEL) || (b_met > 40.)){
    b_step3 = true;
    if (b_step == 2){
      ++b_step;
      if (sys == syst_nom) cutflow_[6][b_channel]++;
    }
  }

  if (selectedJets.size() >1 ){
    b_step4 = true;
    if (b_step == 3){
      ++b_step;
      if (sys == syst_nom) cutflow_[7][b_channel]++;
    }
  }

  if (selectedBJets.size() > 0){
    b_step5 = true;
    if (b_step == 4){
      ++b_step;
      if (sys == syst_nom) cutflow_[8][b_channel]++;
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
      if (jet.bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2M) bjetIndices.push_back(ijet);
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

    if (b_step == 5 and sys == syst_nom) cutflow_[10][b_channel]++;

    if (kinematicReconstructionSolutions.numberOfSolutions()){
      LV top1 = kinematicReconstructionSolutions.solution().top();
      LV top2 = kinematicReconstructionSolutions.solution().antiTop();
	
      b_step7 = true;	
      if (b_step == 5)
	if (sys == syst_nom)
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
    //const bool isBjet1 = jet1->bDiscriminator(BTAG_CSVv2) >= WP_BTAG_CSVv2L;
    const bool isBjet1 = jet1->bDiscriminator(BTAG_CSVv2) >= WP_BTAG_CSVv2M;
    for (auto jet2 = next(jet1); jet2 != end; ++jet2){
      const auto recojet2= jet2->p4();
      //const bool isBjet2 = jet2->bDiscriminator(BTAG_CSVv2) >= WP_BTAG_CSVv2L;
      const bool isBjet2 = jet2->bDiscriminator(BTAG_CSVv2) >= WP_BTAG_CSVv2M;

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
      if (sys == syst_nom) cutflow_[9][b_channel]++;
    }
    //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
    // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_met, b_ll_mass, b_maxweight);
  }
  return true;
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




cat::MuonCollection TtbarDiLeptonAnalyzer::selectMuons(const cat::MuonCollection& muons, systematic sys)
{
  cat::MuonCollection selmuons;
  for (auto& m : muons) {
    cat::Muon mu(m);
    if (std::abs(mu.eta()) > 2.4) continue;
    if (!mu.isTightMuon()) continue;
    
    if (sys == syst_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == syst_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());

    if (mu.pt() < 20.) continue;
    if (mu.relIso(0.4) > 0.15) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }
  return selmuons;
}

cat::ElectronCollection TtbarDiLeptonAnalyzer::selectElecs(const cat::ElectronCollection& elecs, systematic sys)
{
  cat::ElectronCollection selelecs;
  for (auto& e : elecs) {
    cat::Electron el(e);
    if (std::abs(el.eta()) > 2.4) continue;
    if ( !el.isTight() ) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    //if ( !el.electronID("cutBasedElectronID-Summer16-80X-V1-medium") ) continue;
    if ( !el.electronID("cutBasedElectronID-Summer16-80X-V1-tight") ) continue;
    //cout << el.bestTrack()->d0() << endl;
    //if (std::abs(el.bestTrack()->d0()) > 0.0739 ) continue;
    //if (std::abs(el.bestTrack()->dz()) > 0.602 ) continue;
    
    if (sys == syst_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == syst_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if (el.pt() < 20.) continue;
    if (el.relIso(0.3) > 0.0678) continue;
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return selelecs;
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectJets(const cat::JetCollection& jets, const TtbarDiLeptonAnalyzer::LeptonPtrs& recolep, systematic sys)
{
  // Initialize SF_btag
  float Jet_SF_CSV[19];
  for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] = 1.0;

  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    if (sys == syst_jes_u) jet.setP4(j.p4() * j.shiftedEnUp());
    if (sys == syst_jes_d) jet.setP4(j.p4() * j.shiftedEnDown());
    if (sys == syst_jer_u) jet.setP4(j.p4() * j.smearedResUp());
    if (sys == syst_jer_d) jet.setP4(j.p4() * j.smearedResDown());

    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 2.4)  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep->p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    //if (sys == syst_btag_u) b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 1);
    //else if (sys == syst_btag_d) b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, -1);
    //else b_btagweight *= jet.scaleFactorCSVv2(cat::Jet::BTAGCSV_LOOSE, 0);
    for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] *= csvWeight.getSF(jet, iu);
    seljets.push_back(jet);
  }
  for (unsigned int iu=0; iu<19; iu++) b_csvweights.push_back(Jet_SF_CSV[iu]);

  b_btagweight = Jet_SF_CSV[0];
  // if      ( sys == syst_btag_u ) b_btagweight = bTagWeightL.eventWeight(seljets, 1);
  // else if ( sys == syst_btag_d ) b_btagweight = bTagWeightL.eventWeight(seljets, 2);
  // else                          b_btagweight = bTagWeightL.eventWeight(seljets, 0);

  return seljets;
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    //if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2L) continue;
    if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;
    //if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;//forsync
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

void TtbarDiLeptonAnalyzer::setBranch(TTree* tr, systematic sys)
{
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
  
  tr->Branch("tri", &b_tri, "tri/F");
  tr->Branch("tri_up", &b_tri_up, "tri_up/F");
  tr->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  tr->Branch("filtered", &b_filtered, "filtered/O");
  tr->Branch("met", &b_met, "met/F");
  tr->Branch("weight", &b_weight, "weight/F");
  tr->Branch("topPtWeight", &b_topPtWeight, "topPtWeight/F");
  tr->Branch("puweight", &b_puweight, "puweight/F");
  tr->Branch("genweight", &b_genweight, "genweight/F");
  tr->Branch("mueffweight", &b_mueffweight, "mueffweight/F");
  tr->Branch("eleffweight", &b_eleffweight, "eleffweight/F");
  tr->Branch("btagweight", &b_btagweight, "btagweight/F");
  tr->Branch("is3lep", &b_is3lep, "is3lep/I");

  tr->Branch("partonChannel", &b_partonChannel, "partonChannel/I");
  tr->Branch("partonMode1", &b_partonMode1, "partonMode1/I");    
  tr->Branch("partonMode2", &b_partonMode2, "partonMode2/I");    
  tr->Branch("partonInPhase", &b_partonInPhase, "partonInPhase/O");
  tr->Branch("partonInPhaseLep", &b_partonInPhaseLep, "partonInPhaseLep/O");
  tr->Branch("partonInPhaseJet", &b_partonInPhaseJet, "partonInPhaseJet/O");
  tr->Branch("partonlep1_pid", &b_partonlep1_pid, "partonlep1_pid/I");    
  tr->Branch("partonlep2_pid", &b_partonlep2_pid, "partonlep2_pid/I");
  
  tr->Branch("pseudoChannel", &b_pseudoChannel, "pseudoChannel/I");
  tr->Branch("pseudoInPhase", &b_pseudoInPhase, "pseudoInPhase/O");
  tr->Branch("pseudolep1_pid", &b_pseudolep1_pid, "pseudolep1_pid/I");    
  tr->Branch("pseudolep2_pid", &b_pseudolep2_pid, "pseudolep2_pid/I");    

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
  tr->Branch("btagweight_up", &b_btagweight_up, "btagweight_up/F");
  tr->Branch("btagweight_dn", &b_btagweight_dn, "btagweight_dn/F");

  tr->Branch("csvweights","std::vector<float>",&b_csvweights);
  tr->Branch("pdfWeights","std::vector<float>",&b_pdfWeights);
  // tr->Branch("scaleWeights_up","std::vector<float>",&b_scaleWeights_up);
  // tr->Branch("scaleWeights_dn","std::vector<float>",&b_scaleWeights_dn);
  tr->Branch("scaleWeights","std::vector<float>",&b_scaleWeights);

  tr->Branch("partonlep1", "TLorentzVector", &b_partonlep1);
  tr->Branch("partonlep2", "TLorentzVector", &b_partonlep2);
  tr->Branch("partondilep", "TLorentzVector", &b_partondilep);
  tr->Branch("partonjet1", "TLorentzVector", &b_partonjet1);
  tr->Branch("partonjet2", "TLorentzVector", &b_partonjet2);
  tr->Branch("partontop1", "TLorentzVector", &b_partontop1);
  tr->Branch("partontop2", "TLorentzVector", &b_partontop2);
  tr->Branch("partonttbar", "TLorentzVector", &b_partonttbar);
  tr->Branch("partonttbar_dphi", &b_partonttbar_dphi, "partonttbar_dphi/F");

  tr->Branch("pseudolep1", "TLorentzVector", &b_pseudolep1);
  tr->Branch("pseudolep2", "TLorentzVector", &b_pseudolep2);
  tr->Branch("pseudodilep", "TLorentzVector", &b_pseudodilep);
  tr->Branch("pseudojet1", "TLorentzVector", &b_pseudojet1);
  tr->Branch("pseudojet2", "TLorentzVector", &b_pseudojet2);
  tr->Branch("pseudotop1", "TLorentzVector", &b_pseudotop1);
  tr->Branch("pseudotop2", "TLorentzVector", &b_pseudotop2);
  tr->Branch("pseudottbar", "TLorentzVector", &b_pseudottbar);
  tr->Branch("pseudottbar_dphi", &b_pseudottbar_dphi, "pseudottbar_dphi/F");


// ############### Dstar begin #######################

  b_d0         = new TClonesArray("TLorentzVector",100);
  b_d0_dau1    = new TClonesArray("TLorentzVector",100);
  b_d0_dau2    = new TClonesArray("TLorentzVector",100);
  b_dstar      = new TClonesArray("TLorentzVector",100);
  b_dstar_dau1 = new TClonesArray("TLorentzVector",100);
  b_dstar_dau2 = new TClonesArray("TLorentzVector",100);
  b_dstar_dau3 = new TClonesArray("TLorentzVector",100);
  b_Jpsi       = new TClonesArray("TLorentzVector",100);
  b_Jpsi_dau1  = new TClonesArray("TLorentzVector",100);
  b_Jpsi_dau2  = new TClonesArray("TLorentzVector",100);

  // D0
  tr->Branch("d0","TClonesArray",&b_d0,32000,0);
  tr->Branch("d0_dau1","TClonesArray",&b_d0_dau1,32000,0);
  tr->Branch("d0_dau2","TClonesArray",&b_d0_dau2,32000,0);
  tr->Branch("d0_vProb","std::vector<float>",&b_d0_vProb);

  tr->Branch("d0_true","std::vector<bool>",&b_d0_true);
  tr->Branch("d0_fit","std::vector<bool>",&b_d0_fit);

  tr->Branch("d0_L3D","std::vector<float>",&b_d0_L3D);
  tr->Branch("d0_LXY","std::vector<float>",&b_d0_LXY);
  tr->Branch("d0_dRTrue","std::vector<float>",&b_d0_dRTrue);
  tr->Branch("d0_relPtTrue","std::vector<float>",&b_d0_relPtTrue);
  tr->Branch("d0_dca","std::vector<float>",&b_d0_dca);

  tr->Branch("d0_dau1_q","std::vector<float>",&b_d0_dau1_q);
  tr->Branch("d0_dau2_q","std::vector<float>",&b_d0_dau2_q);

  tr->Branch("d0_lepSV_lowM","std::vector<float>",&b_d0_lepSV_lowM);
  tr->Branch("d0_lepSV_dRM","std::vector<float>",&b_d0_lepSV_dRM);
  tr->Branch("d0_lepSV_correctM","std::vector<float>",&b_d0_lepSV_correctM); // for test

  tr->Branch("dstar",    "TClonesArray",&b_dstar    ,32000,0);
  tr->Branch("dstar_dau1","TClonesArray",&b_dstar_dau1,32000,0);
  tr->Branch("dstar_dau2","TClonesArray",&b_dstar_dau2,32000,0);
  tr->Branch("dstar_dau3","TClonesArray",&b_dstar_dau3,32000,0);

  tr->Branch("dstar_true","std::vector<bool>",&b_dstar_true);
  tr->Branch("dstar_fit","std::vector<bool>",&b_dstar_fit);
  tr->Branch("dstar_L3D","std::vector<float>",&b_dstar_L3D);
  tr->Branch("dstar_LXY","std::vector<float>",&b_dstar_LXY);
  tr->Branch("dstar_dRTrue","std::vector<float>",&b_dstar_dRTrue);
  tr->Branch("dstar_relPtTrue","std::vector<float>",&b_dstar_relPtTrue);
  tr->Branch("dstar_dca","std::vector<float>",&b_dstar_dca);
  tr->Branch("dstar_dca2","std::vector<float>",&b_dstar_dca2);
  tr->Branch("dstar_dca3","std::vector<float>",&b_dstar_dca3);

  tr->Branch("dstar_dau1_q","std::vector<float>",&b_dstar_dau1_q);
  tr->Branch("dstar_dau2_q","std::vector<float>",&b_dstar_dau2_q);
  tr->Branch("dstar_dau3_q","std::vector<float>",&b_dstar_dau3_q);
  tr->Branch("dstar_vProb","std::vector<float>",&b_dstar_vProb);
  tr->Branch("dstar_diffMass","std::vector<float>",&b_dstar_diffMass);

  tr->Branch("dstar_lepSV_lowM","std::vector<float>",&b_dstar_lepSV_lowM);
  tr->Branch("dstar_lepSV_dRM","std::vector<float>",&b_dstar_lepSV_dRM);
  tr->Branch("dstar_opCharge_M","std::vector<float>",&b_dstar_opCharge_M);
  tr->Branch("dstar_lepSV_correctM","std::vector<float>",&b_dstar_lepSV_correctM);

  tr->Branch("Jpsi","TClonesArray",&b_Jpsi,32000,0);
  tr->Branch("Jpsi_dau1","TClonesArray",&b_Jpsi_dau1,32000,0);
  tr->Branch("Jpsi_dau2","TClonesArray",&b_Jpsi_dau2,32000,0);
  tr->Branch("Jpsi_vProb","std::vector<float>",&b_Jpsi_vProb);

  tr->Branch("Jpsi_true","std::vector<bool>",&b_Jpsi_true);
  tr->Branch("Jpsi_fit","std::vector<bool>",&b_Jpsi_fit);

  tr->Branch("Jpsi_L3D","std::vector<float>",&b_Jpsi_L3D);
  tr->Branch("Jpsi_LXY","std::vector<float>",&b_Jpsi_LXY);
  tr->Branch("Jpsi_dRTrue","std::vector<float>",&b_Jpsi_dRTrue);
  tr->Branch("Jpsi_relPtTrue","std::vector<float>",&b_Jpsi_relPtTrue);
  tr->Branch("Jpsi_dca","std::vector<float>",&b_Jpsi_dca);

  tr->Branch("Jpsi_dau1_q","std::vector<float>",&b_Jpsi_dau1_q);
  tr->Branch("Jpsi_dau2_q","std::vector<float>",&b_Jpsi_dau2_q);
  tr->Branch("Jpsi_dau_pid","std::vector<int>",&b_Jpsi_dau_pid);

  tr->Branch("Jpsi_lepSV_lowM","std::vector<float>",&b_Jpsi_lepSV_lowM);
  tr->Branch("Jpsi_lepSV_dRM","std::vector<float>",&b_Jpsi_lepSV_dRM);
  tr->Branch("Jpsi_lepSV_correctM","std::vector<float>",&b_Jpsi_lepSV_correctM); // for test


// ############### Dstar end  #######################

}

void TtbarDiLeptonAnalyzer::resetBranch()
{
  b_nvertex = 0;b_step = -1;b_channel = 0;b_njet = 0;b_nbjet = 0;
  b_step1 = 0;b_step2 = 0;b_step3 = 0;b_step4 = 0;b_step5 = 0;b_step6 = 0;b_step7 = 0;b_tri = 0;b_filtered = 0;
  b_met = -9;
  b_weight = 1; b_puweight = 1; b_puweight_up = 1; b_puweight_dn = 1; b_genweight = 1;
  b_mueffweight = 1;b_mueffweight_up = 1;b_mueffweight_dn = 1;
  b_eleffweight = 1;b_eleffweight_up = 1;b_eleffweight_dn = 1;
  b_btagweight = 1;b_btagweight_up = 1;b_btagweight_dn = 1;
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

  b_d0->Clear();    b_d0_dau1->Clear();    b_d0_dau2->Clear();
  b_dstar->Clear(); b_dstar_dau1->Clear(); b_dstar_dau2->Clear(); b_dstar_dau3->Clear();
  b_Jpsi->Clear();    b_Jpsi_dau1->Clear();    b_Jpsi_dau2->Clear();

// ##################### Dstar begin ######################
  // D0
  b_d0_true.clear() ;
  b_d0_LXY.clear(); b_d0_L3D.clear(); b_d0_fit.clear(); b_d0_dRTrue.clear(); b_d0_relPtTrue.clear(); b_d0_dca.clear();

  b_d0_dau1_q.clear();
  b_d0_dau2_q.clear();
  b_d0_vProb.clear();

  b_d0_lepSV_lowM.clear();
  b_d0_lepSV_dRM.clear();
  b_d0_lepSV_correctM.clear();

  // Dstar
  b_dstar_true.clear();
  b_dstar_LXY.clear(); b_dstar_L3D.clear(); b_dstar_fit.clear(); b_dstar_dRTrue.clear(); b_dstar_relPtTrue.clear(); b_dstar_dca.clear(); b_dstar_dca2.clear(); b_dstar_dca3.clear();

  b_dstar_dau1_q.clear();
  b_dstar_dau2_q.clear();
  b_dstar_dau3_q.clear();
  b_dstar_vProb.clear();
  b_dstar_diffMass.clear();

  b_dstar_lepSV_lowM.clear();
  b_dstar_lepSV_dRM.clear();
  b_dstar_opCharge_M.clear();
  b_dstar_lepSV_correctM.clear();
  // Jpsi
  b_Jpsi_true.clear() ;
  b_Jpsi_LXY.clear(); b_Jpsi_L3D.clear(); b_Jpsi_fit.clear(); b_Jpsi_dRTrue.clear(); b_Jpsi_relPtTrue.clear(); b_Jpsi_dca.clear();

  b_Jpsi_dau1_q.clear();
  b_Jpsi_dau2_q.clear();
  b_Jpsi_dau_pid.clear();
  b_Jpsi_vProb.clear();

  b_Jpsi_lepSV_lowM.clear();
  b_Jpsi_lepSV_dRM.clear();
  b_Jpsi_lepSV_correctM.clear();


// ############################ Dstar end ##########################
}

bool TtbarDiLeptonAnalyzer::genInformation(const edm::Event& iEvent)
{
  bool keepTtbarSignal = false;

  edm::Handle<int> partonTop_channel;
  if ( iEvent.getByToken(partonTop_channel_, partonTop_channel)){

    edm::Handle<float> topPtWeightHandle;
    iEvent.getByToken(topPtWeight_, topPtWeightHandle);
    b_topPtWeight = *topPtWeightHandle;

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
  return keepTtbarSignal;
}
//define this as a plug-in
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
