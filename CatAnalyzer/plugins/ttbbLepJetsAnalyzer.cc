// -*- C++ -*-
//
// Package:    ttbbLepJetsAnalyzer
// Class:      ttbbLepJetsAnalyzer
//
/**\class ttbbLepJetsAnalyzer ttbbLepJetsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Javier Brochero Cifuentes,512 1-001,+41789428316
//         Created:  Tue Feb  3 09:52:55 CET 2015
// $Id$
//
//

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenWeights.h"

#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Kinematic Reconstruction
#include "CATTools/CatAnalyzer/interface/LepJets_Fitter.h"

#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//


using namespace cat;

class ttbbLepJetsAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ttbbLepJetsAnalyzer(const edm::ParameterSet&);
  ~ttbbLepJetsAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  //----------------------------------------------------------------
  bool IsSelectMuon    (const cat::Muon     & i_muon_candidate);
  bool IsVetoMuon      (const cat::Muon     & i_muon_candidate);
  bool IsSelectElectron(const cat::Electron & i_electron_candidate);
  bool IsVetoElectron  (const cat::Electron & i_electron_candidate);

  bool isMC_ ;

  int TTbarMC_; // 0->No ttbar, 1->ttbar Signal, 2->ttbar Background
  int TTbarCatMC_;
  unsigned int SkimNJets_;
  bool KFUsebtag_;
  bool CSVPosConKF_;
  // Trigger Names
  std::vector<string> triggerNameDataEl_;
  std::vector<string> triggerNameDataMu_;
  std::vector<string> triggerNameMCEl_;
  std::vector<string> triggerNameMCMu_;
  
  // Event Weights
  edm::EDGetTokenT<float>                        genWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           pdfWeightsToken_;
  edm::EDGetTokenT<std::vector<float>>           scaleUpWeightsToken_;
  edm::EDGetTokenT<std::vector<float>>           scaleDownWeightsToken_;
  edm::EDGetTokenT<float>                        puWeightToken_;
  edm::EDGetTokenT<float>                        puUpWeightToken_;
  edm::EDGetTokenT<float>                        puDownWeightToken_;
  // Object Collections
  edm::EDGetTokenT<reco::GenParticleCollection>  genToken_;
  edm::EDGetTokenT<reco::GenJetCollection>       genJetToken_;
  edm::EDGetTokenT<int>                          genttbarHiggsCatToken_;
  edm::EDGetTokenT<cat::GenTopCollection>        genttbarCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
  edm::EDGetTokenT<cat::JetCollection>           jetToken_;
  edm::EDGetTokenT<cat::METCollection>           metToken_;
  edm::EDGetTokenT<int>                          pvToken_;
  edm::EDGetTokenT<std::vector<vector<int>>>     JetMotherToken_;
  // Trigger
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits2_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

// ----------member data ---------------------------

  TTree *tree;
  TTree *gentree;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Tree Branches
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Event info
  int b_Event, b_Run, b_Lumi_Number;
  float b_GenWeight;
  std::vector<float> *b_ScaleWeight;
  // PU/Vertices
  std::vector<float> *b_PUWeight;
  int b_nGoodPV;
  // Channel and Categorization
  int b_GenChannel;
  int b_Channel;
  std::vector<int>   *b_GenConeCatID;
  std::vector<float> *b_GenCone_gJet_pT;
  std::vector<float> *b_GenCone_gJet_eta;
  std::vector<float> *b_GenCone_gJet_phi;
  std::vector<float> *b_GenCone_gJet_E;
  std::vector<int>   *b_GenCone_gJetFlavW;
  std::vector<int>   *b_GenCone_gJetIndex;
  int b_GenCone_NgJetsW;

  int b_GenHiggsCatID;
  // MET
  float b_MET, b_MET_phi;
  // GEN Leptons
  float b_GenLepton_pT;
  float b_GenLepton_eta;
  float b_GenLepton_phi;
  float b_GenLepton_E;
  // GEN Neutrino
  float b_GenNu_pT;
  float b_GenNu_eta;
  float b_GenNu_phi;
  float b_GenNu_E;

  // Leptons
  float b_Lepton_pT;
  float b_Lepton_eta;
  float b_Lepton_phi;
  float b_Lepton_E;
  std::vector<float> *b_Lepton_SF;
  float b_Lepton_LES;
  // GEN Jets
  std::vector<float> *b_GenJet_pT;
  std::vector<float> *b_GenJet_eta;
  std::vector<float> *b_GenJet_phi;
  std::vector<float> *b_GenJet_E;
  // Jet Mother (MC Studies)
  std::vector<int> *b_GenJet_mom, *b_GenJet_GenConeMom;
  // Jets
  int b_Jet_Number;
  std::vector<float> *b_Jet_pT;
  std::vector<float> *b_Jet_eta;
  std::vector<float> *b_Jet_phi;
  std::vector<float> *b_Jet_E;
  std::vector<int>   *b_Jet_Index, *b_Jet_GenConeMom;
  std::vector<int>   *b_Jet_MatchedGenJetIndex;
  float               b_DRAddJets;
  // Jet Flavour
  std::vector<int> *b_Jet_partonFlavour;
  std::vector<int> *b_Jet_hadronFlavour;
  // JES and JER
  std::vector<float> *b_Jet_JES_Up;
  std::vector<float> *b_Jet_JES_Down;
  std::vector<float> *b_Jet_JER_Up;
  std::vector<float> *b_Jet_JER_Nom;
  std::vector<float> *b_Jet_JER_Down;
  // b-Jet discriminant
  std::vector<float> *b_Jet_CSV;
  std::vector<float> *b_Jet_SF_CSV_25;
  std::vector<float> *b_Jet_SF_CSV_30;
  std::vector<float> *b_Jet_SF_CSV_35;
  std::vector<float> *b_Jet_SF_CSV_40;
  // c-Jet discriminant
  std::vector<float> *b_Jet_CvsL, *b_Jet_CvsB;

  // Kinematic Reconstruction
  float b_Kin_Chi2;

  // KR Leptons
  float b_KinNu_pT;
  float b_KinNu_eta;
  float b_KinNu_phi;
  float b_KinNu_E;
  // KR Jets
  std::vector<float> *b_KinJet_pT;
  std::vector<float> *b_KinJet_eta;
  std::vector<float> *b_KinJet_phi;
  std::vector<float> *b_KinJet_E;

  std::vector<int>   *b_KinJet_Index;


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Histograms
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Histograms: Number of Events and Weights
  TH1F *EventInfo, *ScaleWeights;
  // Scale factor evaluators
  BTagWeightEvaluator SF_CSV_;
  ScaleFactorEvaluator SF_muon_, SF_elec_;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ttbbLepJetsAnalyzer::ttbbLepJetsAnalyzer(const edm::ParameterSet& iConfig):
  TTbarMC_    (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0)),
  TTbarCatMC_ (iConfig.getUntrackedParameter<int>("TTbarCatLabel", 0)),
  SkimNJets_  (iConfig.getUntrackedParameter<unsigned int>("Skim_N_Jets", 0)),
  KFUsebtag_  (iConfig.getUntrackedParameter<bool>("KFUsebtagLabel", true)),
  CSVPosConKF_(iConfig.getUntrackedParameter<bool>("CSVPosConKFLabel", true)),
  triggerNameDataEl_(iConfig.getUntrackedParameter<std::vector<string>>("triggerNameDataEl")),
  triggerNameDataMu_(iConfig.getUntrackedParameter<std::vector<string>>("triggerNameDataMu")),
  triggerNameMCEl_  (iConfig.getUntrackedParameter<std::vector<string>>("triggerNameMCEl")),
  triggerNameMCMu_  (iConfig.getUntrackedParameter<std::vector<string>>("triggerNameMCMu"))
{
  const auto elecSFSet = iConfig.getParameter<edm::ParameterSet>("elecSF");
  SF_elec_.set(elecSFSet.getParameter<std::vector<double>>("pt_bins" ),
	       elecSFSet.getParameter<std::vector<double>>("abseta_bins"),
	       elecSFSet.getParameter<std::vector<double>>("values"  ),
	       elecSFSet.getParameter<std::vector<double>>("errors"  ));
  
  const auto muonSFSet = iConfig.getParameter<edm::ParameterSet>("muonSF");
  SF_muon_.set(muonSFSet.getParameter<std::vector<double>>("pt_bins"    ),
	       muonSFSet.getParameter<std::vector<double>>("abseta_bins"),
	       muonSFSet.getParameter<std::vector<double>>("values"     ),
	       muonSFSet.getParameter<std::vector<double>>("errors"     ));
  
  SF_CSV_.initCSVWeight(false, "csvv2");
  
  // Weights
  genWeightToken_        = consumes<float>             (iConfig.getParameter<edm::InputTag>("genWeightLabel"));
  pdfWeightsToken_       = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("pdfWeightLabel"));
  scaleUpWeightsToken_   = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("scaleUpWeightLabel"));
  scaleDownWeightsToken_ = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("scaleDownWeightLabel"));
  puWeightToken_         = consumes<float>             (iConfig.getParameter<edm::InputTag>("puWeightLabel"));
  puUpWeightToken_       = consumes<float>             (iConfig.getParameter<edm::InputTag>("puUpWeightLabel"));
  puDownWeightToken_     = consumes<float>             (iConfig.getParameter<edm::InputTag>("puDownWeightLabel"));
  // GEN and ttbar Categorization
  genToken_              = consumes<reco::GenParticleCollection>  (iConfig.getParameter<edm::InputTag>("genLabel"));
  genJetToken_           = consumes<reco::GenJetCollection>       (iConfig.getParameter<edm::InputTag>("genJetLabel"));  
  genttbarCatToken_      = consumes<cat::GenTopCollection>        (iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));
  genttbarHiggsCatToken_ = consumes<int>                          (iConfig.getParameter<edm::InputTag>("genHiggsCatLabel"));
  // Object Collections
  muonToken_         = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_     = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_          = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_          = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_           = consumes<int>                          (iConfig.getParameter<edm::InputTag>("pvLabel"));
  JetMotherToken_    = consumes<vector<vector<int>>>          (iConfig.getParameter<edm::InputTag>("JetMother"));
  // Trigger  
  triggerBits_       = consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerBits2_      = consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerBits2"));
  triggerObjects_    = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerObjects"));
  
  b_PUWeight     = new std::vector<float>;
  b_ScaleWeight  = new std::vector<float>;  
  b_Lepton_SF    = new std::vector<float>;  

  b_GenConeCatID      = new std::vector<int>;
  b_GenCone_gJet_pT   = new std::vector<float>;
  b_GenCone_gJet_eta  = new std::vector<float>;
  b_GenCone_gJet_phi  = new std::vector<float>;
  b_GenCone_gJet_E    = new std::vector<float>;
  b_GenCone_gJetIndex = new std::vector<int>;
  b_GenCone_gJetFlavW = new std::vector<int>;

  b_GenJet_pT = new std::vector<float>;
  b_GenJet_eta = new std::vector<float>;
  b_GenJet_phi = new std::vector<float>;
  b_GenJet_E  = new std::vector<float>;
  b_GenJet_mom = new std::vector<int>;
  b_GenJet_GenConeMom = new std::vector<int>;
  b_Jet_MatchedGenJetIndex = new std::vector<int>;
  b_Jet_pT   = new std::vector<float>;
  b_Jet_eta   = new std::vector<float>;
  b_Jet_phi   = new std::vector<float>;
  b_Jet_E    = new std::vector<float>;
  b_Jet_Index= new std::vector<int>;
  b_Jet_GenConeMom=new std::vector<int>;
  b_Jet_CSV  = new std::vector<float>;
  b_Jet_SF_CSV_25  = new std::vector<float>;
  b_Jet_SF_CSV_30  = new std::vector<float>;
  b_Jet_SF_CSV_35  = new std::vector<float>;
  b_Jet_SF_CSV_40  = new std::vector<float>;
  b_Jet_CvsL    = new std::vector<float>;  
  b_Jet_CvsB    = new std::vector<float>;  
  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;
  b_Jet_JES_Up   = new std::vector<float>;
  b_Jet_JES_Down = new std::vector<float>;
  b_Jet_JER_Up   = new std::vector<float>;
  b_Jet_JER_Nom  = new std::vector<float>;
  b_Jet_JER_Down = new std::vector<float>;
  // Kinematic Reconstruction
  b_KinJet_pT  = new std::vector<float>;
  b_KinJet_eta = new std::vector<float>;
  b_KinJet_phi = new std::vector<float>;
  b_KinJet_E   = new std::vector<float>;

  b_KinJet_Index = new std::vector<int>;


  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "Event/I");
  tree->Branch("run",        &b_Run,         "Run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "Lumi_Number/I");
  tree->Branch("genweight",  &b_GenWeight,   "GenWeight/F");
  tree->Branch("GoodPV",     &b_nGoodPV,     "nGoodPV/I");
  tree->Branch("channel",    &b_Channel,     "Channel/I");

  tree->Branch("PUWeight",   "std::vector<float>", &b_PUWeight);

  tree->Branch("MET",     &b_MET,     "MET/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("lepton_pT",  &b_Lepton_pT,  "lepton_pT/F");
  tree->Branch("lepton_eta", &b_Lepton_eta,  "lepton_eta/F");
  tree->Branch("lepton_phi", &b_Lepton_phi, "lepton_phi/F");
  tree->Branch("lepton_E" ,  &b_Lepton_E,   "lepton_E/F" );
  tree->Branch("lepton_LES", &b_Lepton_LES, "lepton_LES/F" );

  tree->Branch("lepton_SF",  "std::vector<float>", &b_Lepton_SF );

  tree->Branch("jet_pT",           "std::vector<float>", &b_Jet_pT);
  tree->Branch("jet_eta",          "std::vector<float>", &b_Jet_eta);
  tree->Branch("jet_phi",          "std::vector<float>", &b_Jet_phi);
  tree->Branch("jet_E" ,           "std::vector<float>", &b_Jet_E );
  tree->Branch("jet_index" ,       "std::vector<int>",   &b_Jet_Index );
  tree->Branch("jet_gencone_mom" , "std::vector<int>",   &b_Jet_GenConeMom );
  tree->Branch("jet_CSV" ,         "std::vector<float>", &b_Jet_CSV );
  tree->Branch("jet_SF_CSV_25",    "std::vector<float>", &b_Jet_SF_CSV_25 );
  tree->Branch("jet_SF_CSV_30",    "std::vector<float>", &b_Jet_SF_CSV_30 );
  tree->Branch("jet_SF_CSV_35",    "std::vector<float>", &b_Jet_SF_CSV_35 );
  tree->Branch("jet_SF_CSV_40",    "std::vector<float>", &b_Jet_SF_CSV_40 );
  tree->Branch("jet_CvsL",         "std::vector<float>", &b_Jet_CvsL );
  tree->Branch("jet_CvsB",         "std::vector<float>", &b_Jet_CvsB );

  tree->Branch("jet_Number" , &b_Jet_Number, "jet_number/I" );

  tree->Branch("jet_partonFlavour", "std::vector<int>",   &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>",   &b_Jet_hadronFlavour);
  tree->Branch("jet_JES_Up",        "std::vector<float>", &b_Jet_JES_Up );
  tree->Branch("jet_JES_Down",      "std::vector<float>", &b_Jet_JES_Down );
  tree->Branch("jet_JER_Up",        "std::vector<float>", &b_Jet_JER_Up );
  tree->Branch("jet_JER_Nom",       "std::vector<float>", &b_Jet_JER_Nom );
  tree->Branch("jet_JER_Down",      "std::vector<float>", &b_Jet_JER_Down );

  // Kinematic Reconstruction
  tree->Branch("kin_chi2",   &b_Kin_Chi2,  "kin_chi2/F");
  tree->Branch("kinnu_pT",   &b_KinNu_pT,  "kinnu_pT/F");
  tree->Branch("kinnu_eta",  &b_KinNu_eta, "kinnu_eta/F");
  tree->Branch("kinnu_phi",  &b_KinNu_phi, "kinnu_phi/F");
  tree->Branch("kinnu_E",    &b_KinNu_E,   "kinnu_E/F");
  
  tree->Branch("kinjet_pT",    "std::vector<float>", &b_KinJet_pT);
  tree->Branch("kinjet_eta",   "std::vector<float>", &b_KinJet_eta);
  tree->Branch("kinjet_phi",   "std::vector<float>", &b_KinJet_phi);
  tree->Branch("kinjet_E",     "std::vector<float>", &b_KinJet_E);
  tree->Branch("kinjet_index", "std::vector<int>",   &b_KinJet_Index);
  

  // GEN Variables (only ttbarSignal)
  if(TTbarMC_ == 1){
    tree->Branch("scaleweight",   "std::vector<float>", &b_ScaleWeight );

    tree->Branch("jet_MatchedGenJetIndex", "std::vector<int>",  &b_Jet_MatchedGenJetIndex);

    tree->Branch("genconecatid" ,      "std::vector<int>",   &b_GenConeCatID);
    tree->Branch("gencone_gjet_pT" ,   "std::vector<float>", &b_GenCone_gJet_pT);
    tree->Branch("gencone_gjet_eta" ,  "std::vector<float>", &b_GenCone_gJet_eta);
    tree->Branch("gencone_gjet_phi" ,  "std::vector<float>", &b_GenCone_gJet_phi);
    tree->Branch("gencone_gjet_E" ,    "std::vector<float>", &b_GenCone_gJet_E);
    tree->Branch("gencone_gjetIndex" , "std::vector<int>",   &b_GenCone_gJetIndex);
    tree->Branch("gencone_gJetFlavW" , "std::vector<int>",   &b_GenCone_gJetFlavW);

    tree->Branch("gencone_NgjetsW", &b_GenCone_NgJetsW, "b_GenCone_NgJetsW/I");
    tree->Branch("draddjets",       &b_DRAddJets,       "DRAddJets/F");
    tree->Branch("genhiggscatid",   &b_GenHiggsCatID,   "genhiggscatid/I");

    tree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    tree->Branch("genlepton_pT",  &b_GenLepton_pT,  "genlepton_pT/F");
    tree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    tree->Branch("genlepton_phi", &b_GenLepton_phi, "genlepton_phi/F");
    tree->Branch("genlepton_E",   &b_GenLepton_E,   "genlepton_E/F");
    tree->Branch("gennu_pT",      &b_GenNu_pT,      "gennu_pT/F");
    tree->Branch("gennu_eta",     &b_GenNu_eta,     "gennu_eta/F");
    tree->Branch("gennu_phi",     &b_GenNu_phi,     "gennu_phi/F");
    tree->Branch("gennu_E",       &b_GenNu_E,       "gennu_E/F");

    tree->Branch("genjet_pT",  "std::vector<float>", &b_GenJet_pT);
    tree->Branch("genjet_eta", "std::vector<float>", &b_GenJet_eta);
    tree->Branch("genjet_phi", "std::vector<float>", &b_GenJet_phi);
    tree->Branch("genjet_E",   "std::vector<float>", &b_GenJet_E);
    tree->Branch("genjet_mom", "std::vector<int>",  &b_GenJet_mom);

    tree->Branch("genjet_gencone_mom", "std::vector<int>",  &b_GenJet_GenConeMom);

    //GEN TREE
    gentree = fs->make<TTree>("gentree", "TopGENTree");
    gentree->Branch("genweight",     &b_GenWeight,     "GenWeight/F");
    gentree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    gentree->Branch("genhiggscatid", &b_GenHiggsCatID, "genhiggscatid/I");
    gentree->Branch("draddjets",     &b_DRAddJets,     "DRAddJets/F");
    gentree->Branch("genlepton_pT",  &b_GenLepton_pT,  "genlepton_pT/F");
    gentree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");

    gentree->Branch("scaleweight",        "std::vector<float>", &b_ScaleWeight );
    gentree->Branch("genconecatid" ,      "std::vector<int>",   &b_GenConeCatID);
    gentree->Branch("genjet_pT",          "std::vector<float>", &b_GenJet_pT);
    gentree->Branch("genjet_eta",         "std::vector<float>", &b_GenJet_eta);
    gentree->Branch("genjet_phi",         "std::vector<float>", &b_GenJet_phi);
    gentree->Branch("genjet_E",           "std::vector<float>", &b_GenJet_E);
    gentree->Branch("genjet_mom",         "std::vector<int>",   &b_GenJet_mom);
    gentree->Branch("genjet_gencone_mom", "std::vector<int>",   &b_GenJet_GenConeMom);
  }

  EventInfo = fs->make<TH1F>("EventInfo","Event Information",8,0,8);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
  EventInfo->GetXaxis()->SetBinLabel(3,"ttbb Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(4,"ttbb Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(5,"ttbb Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(6,"ttjj Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(7,"ttjj Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(8,"ttjj Lep (tau lep decay) Events pT>20");

  ScaleWeights = fs->make<TH1F>("ScaleWeights","Event Weights:",6,0,6);
  ScaleWeights->GetXaxis()->SetBinLabel(1,"muR=Nom  muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(2,"muR=Nom  muF=Down");
  ScaleWeights->GetXaxis()->SetBinLabel(3,"muR=Up   muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(4,"muR=Up   muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(5,"muR=Down muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(6,"muR=Down muF=Down");


}


ttbbLepJetsAnalyzer::~ttbbLepJetsAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete b_PUWeight;
  delete b_ScaleWeight;

  delete b_GenConeCatID;
  delete b_GenCone_gJet_pT;
  delete b_GenCone_gJet_eta;
  delete b_GenCone_gJet_phi;
  delete b_GenCone_gJet_E;
  delete b_GenCone_gJetIndex;
  delete b_GenCone_gJetFlavW;

  delete b_GenJet_pT;
  delete b_GenJet_eta;
  delete b_GenJet_phi;
  delete b_GenJet_E;
  delete b_GenJet_mom;
  delete b_GenJet_GenConeMom;
  delete b_Lepton_SF;

  delete b_Jet_pT;
  delete b_Jet_eta;
  delete b_Jet_phi;
  delete b_Jet_E;
  delete b_Jet_Index;
  delete b_Jet_GenConeMom;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_MatchedGenJetIndex;

  delete b_Jet_JES_Up;
  delete b_Jet_JES_Down;
  delete b_Jet_JER_Up;
  delete b_Jet_JER_Nom;
  delete b_Jet_JER_Down;

  delete b_Jet_SF_CSV_25;
  delete b_Jet_SF_CSV_30;
  delete b_Jet_SF_CSV_35;
  delete b_Jet_SF_CSV_40;

  delete b_Jet_CSV;
  delete b_Jet_CvsL;
  delete b_Jet_CvsB;

  delete b_KinJet_pT;
  delete b_KinJet_eta;
  delete b_KinJet_phi;
  delete b_KinJet_E;
  
  delete b_KinJet_Index;

}

//
// member functions
//

// ------------ method called for each event  ------------
void ttbbLepJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_PUWeight   ->clear();
  b_ScaleWeight->clear();

  b_GenConeCatID->clear();
  b_GenCone_gJet_pT->clear();
  b_GenCone_gJet_eta->clear();
  b_GenCone_gJet_phi->clear();
  b_GenCone_gJet_E->clear();
  b_GenCone_gJetIndex->clear();
  b_GenCone_gJetFlavW->clear();


  b_GenJet_pT  ->clear();
  b_GenJet_eta ->clear();
  b_GenJet_phi ->clear();
  b_GenJet_E   ->clear();
  b_GenJet_mom->clear();
  b_GenJet_GenConeMom->clear();

  b_Lepton_SF->clear();

  b_Jet_pT    ->clear();
  b_Jet_eta   ->clear();
  b_Jet_phi   ->clear();
  b_Jet_E     ->clear();
  b_Jet_Index->clear();
  b_Jet_GenConeMom->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_MatchedGenJetIndex->clear();

  b_Jet_JES_Up  ->clear();
  b_Jet_JES_Down->clear();
  b_Jet_JER_Up  ->clear();
  b_Jet_JER_Nom ->clear();
  b_Jet_JER_Down->clear();

  b_Jet_CSV      ->clear();
  b_Jet_SF_CSV_25->clear();
  b_Jet_SF_CSV_30->clear();
  b_Jet_SF_CSV_35->clear();
  b_Jet_SF_CSV_40->clear();
  b_Jet_CvsL     ->clear();
  b_Jet_CvsB     ->clear();

  b_KinJet_pT ->clear();
  b_KinJet_eta->clear();
  b_KinJet_phi->clear();
  b_KinJet_E  ->clear();
  
  b_KinJet_Index->clear();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  isMC_ = !iEvent.isRealData();

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();

  if(isMC_) {

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // PU Info
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    edm::Handle<float> PUWeight;
    iEvent.getByToken(puWeightToken_, PUWeight);
    b_PUWeight->push_back(*PUWeight); // Central

    edm::Handle<float> PUWeight_Up;
    iEvent.getByToken(puUpWeightToken_, PUWeight_Up);
    b_PUWeight->push_back(*PUWeight_Up); //Syst. Up

    edm::Handle<float> PUWeight_Down;
    iEvent.getByToken(puDownWeightToken_, PUWeight_Down);
    b_PUWeight->push_back(*PUWeight_Down); //Syst. Up

    //---------------------------------------------------------------------------
    // Weights at Generation Level: MC@NLO
    //---------------------------------------------------------------------------

    edm::Handle<float> genWeight;

    iEvent.getByToken(genWeightToken_, genWeight);
    b_GenWeight = *genWeight;

    EventInfo->Fill(0.5, 1.0);         // Number of Events
    EventInfo->Fill(1.5, b_GenWeight); // Sum of Weights
  }

  else{
    b_PUWeight->push_back(1.0);
    b_GenWeight = 1.0;
  }

  //---------------------------------------------------------------------------
  // Weights for Syst. Scale: ttbar
  //---------------------------------------------------------------------------
  if(TTbarMC_ == 1 ) {
    edm::Handle<std::vector<float>> scaleUpWeightsHandle, scaleDownWeightsHandle;
    iEvent.getByToken(scaleUpWeightsToken_,   scaleUpWeightsHandle);
    iEvent.getByToken(scaleDownWeightsToken_, scaleDownWeightsHandle);

    // muR/muF Scale Weights
    b_ScaleWeight->push_back(scaleUpWeightsHandle  ->at(0)); // muR=Nom  muF=Up
    b_ScaleWeight->push_back(scaleDownWeightsHandle->at(0)); // muR=Nom  muF=Down
    b_ScaleWeight->push_back(scaleUpWeightsHandle  ->at(1)); // muR=Up   muF=Nom
    b_ScaleWeight->push_back(scaleUpWeightsHandle  ->at(2)); // muR=Up   muF=Up
    b_ScaleWeight->push_back(scaleDownWeightsHandle->at(1)); // muR=Down muF=Nom
    b_ScaleWeight->push_back(scaleDownWeightsHandle->at(2)); // muR=Down muF=Down
    
    // Sum of muR/muF Scale Weights
    for(unsigned int iscale = 0; iscale< b_ScaleWeight->size(); iscale++)
      ScaleWeights->Fill(iscale, b_ScaleWeight->at(iscale)); 
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Generated Particles (For Pythia8)
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  int nGenLep  = 999;

  bool IsCat = false;

  if(TTbarMC_ > 0) {

    //---------------------------------------------------------------------------
    // Event Categorization Using Higgs Code
    // // Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GenHFHadronMatcher
    //---------------------------------------------------------------------------
    edm::Handle<int> genttbarHiggsCatHandle;
    iEvent.getByToken( genttbarHiggsCatToken_, genttbarHiggsCatHandle );
    b_GenHiggsCatID = *genttbarHiggsCatHandle;

    //---------------------------------------------------------------------------
    // Event Categorization Using Cone
    //---------------------------------------------------------------------------
    edm::Handle<cat::GenTopCollection> genttbarConeCat;
    iEvent.getByToken( genttbarCatToken_, genttbarConeCat );

    // [0]: Decay mode. Just as a Cross Check!
    b_GenConeCatID->push_back(genttbarConeCat->begin()->semiLeptonic(0));
    // [1]: Number of Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NJets20());
    // [2]: Number of b-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NbJets20());
    // [3]: Number of c-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NcJets20());
    // [4]: Number of b-Jets Not comming from the top
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NbJets20NoTop());
    // [5]: Number of add Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()-> NaddJets20());
    // [6]: Number of add b-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()-> NaddbJets20());
    // [7]: Number of add c-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()-> NaddcJets20());

    // additional Jets Information
    // p4
    math::XYZTLorentzVector gJetGenCone[6];
    gJetGenCone[0] = genttbarConeCat->begin()->bJetsFromTop1();
    gJetGenCone[1] = genttbarConeCat->begin()->bJetsFromTop2();
    gJetGenCone[2] = genttbarConeCat->begin()->JetsFromW1();
    gJetGenCone[3] = genttbarConeCat->begin()->JetsFromW2();
    gJetGenCone[4] = genttbarConeCat->begin()->addJets1();
    gJetGenCone[5] = genttbarConeCat->begin()->addJets2();

    b_GenCone_NgJetsW = genttbarConeCat->begin()-> NWJets();
    
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()-> Wquarkflav1());
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()-> Wquarkflav2());

    for (int ijGT=0; ijGT<6; ijGT++){
      b_GenCone_gJet_pT  ->push_back(gJetGenCone[ijGT].Pt());
      b_GenCone_gJet_eta ->push_back(gJetGenCone[ijGT].Eta());
      b_GenCone_gJet_phi ->push_back(gJetGenCone[ijGT].Phi());
      b_GenCone_gJet_E   ->push_back(gJetGenCone[ijGT].E());
    }

    // DR 
    b_DRAddJets = genttbarConeCat->begin()->dRaddJets();

    if(genttbarConeCat->begin()-> NaddbJets20() > 1){
      EventInfo->Fill(2.5, 1.0); // Number of ttbb Events	
      if(genttbarConeCat->begin()->semiLeptonic(-1)) EventInfo->Fill(3.5, 1.0); // Number of ttbb Events (Includes tau) 
      if(genttbarConeCat->begin()->semiLeptonic(0))  EventInfo->Fill(4.5, 1.0); // Number of ttbb Events (Includes tau leptonic decay) 
    }
    if(genttbarConeCat->begin()-> NaddJets20() > 1){
      EventInfo->Fill(5.5, 1.0); // Number of ttjj Events	
      if(genttbarConeCat->begin()->semiLeptonic(-1)) EventInfo->Fill(6.5, 1.0); // Number of ttjj Events (Includes tau) 
      if(genttbarConeCat->begin()->semiLeptonic(0))  EventInfo->Fill(7.5, 1.0); // Number of ttjj Events (Includes tau leptonic decay) 
    }
    
    
    //---------------------------------------------------------------------------
    // Using the GenChannel from GenTop categorization
    //---------------------------------------------------------------------------
    nGenLep = genttbarConeCat->begin()->semiLeptonic(0); // semiLeptonic(0) includes tau leptonic decay
    //---------------------------------------------------------------------------

    // Category
    bool Isttjj = false;
    bool Isttbb = false;
    bool Isttb  = false;
    bool Isttcc = false;
    bool IsttLF = false;
    bool Istt   = false;

    // Categorization based in the Full Ph-Sp
    if(genttbarConeCat->begin()->NaddJets20() > 1) Isttjj = true;

    if      (genttbarConeCat->begin()->NaddbJets20() > 1) Isttbb = true;
    else if (genttbarConeCat->begin()->NaddbJets20() > 0) Isttb  = true;
    else if (genttbarConeCat->begin()->NaddcJets20() > 1) Isttcc = true;
    else if (genttbarConeCat->begin()->NaddJets20()  > 1) IsttLF = true;
    else Istt = true;

    // Categorization based in the Visible Ph-Sp
    // if(genttbarConeCat->begin()->NbJets20() > 1 && 
    //    genttbarConeCat->begin()->NJets20()  > 5) Isttjj = true;

    // if      (genttbarConeCat->begin()->NbJets20() > 3  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5) Isttbb = true;
    // else if (genttbarConeCat->begin()->NbJets20() > 2  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5) Isttb  = true;
    // else if (genttbarConeCat->begin()->NbJets20() > 1  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5  && 
    // 	     genttbarConeCat->begin()->NcJets20() > 1) Isttcc = true;
    // else if (genttbarConeCat->begin()->NbJets20() > 1  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5) IsttLF = true;
    // else Istt = true;


    if (TTbarCatMC_ == 0)           IsCat = true;
    if (Isttbb && TTbarCatMC_ == 1) IsCat = true;
    if (Isttb  && TTbarCatMC_ == 2) IsCat = true;
    if (Isttcc && TTbarCatMC_ == 3) IsCat = true;
    if (IsttLF && TTbarCatMC_ == 4) IsCat = true;
    if (Istt   && TTbarCatMC_ == 5) IsCat = true;
    if (Isttjj && TTbarCatMC_ == 6) IsCat = true;


    if(isMC_ && TTbarMC_== 1 && IsCat){
      if(nGenLep == 1){

	if (genttbarConeCat->begin()->semiLeptonicMuo())      b_GenChannel = 0;
	else if (genttbarConeCat->begin()->semiLeptonicEle()) b_GenChannel = 1;
	
	if(genttbarConeCat->begin()->lepton1().pt() != 0.){
	  b_GenLepton_pT  = genttbarConeCat->begin()->lepton1().pt();
	  b_GenLepton_eta = genttbarConeCat->begin()->lepton1().eta();
	}
	else{
	  b_GenLepton_pT  = genttbarConeCat->begin()->lepton2().pt();
	  b_GenLepton_eta = genttbarConeCat->begin()->lepton2().eta();
	}

	//---------------------------------------------------------------------------
	// Generated particles: Lepton, neutrino and jets
	//---------------------------------------------------------------------------
	
	edm::Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByToken(genToken_, genParticles);
	
	TLorentzVector Genlepton(0,0,0,0); 
	TLorentzVector Gennu(0,0,0,0);

	// Gen Status: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
	// abs(id) == 6  // t-Quark
	// abs(id) == 5  // b-Quark
	// abs(id) == 24 // W-Boson
	// abs(id) == 15 // Tau
	// abs(id) == 11 // Electron
	// abs(id) == 13 // Muon

	for (unsigned int i = 0; i < genParticles->size(); i++){

	  const reco::Candidate & gp = (*genParticles)[i];
	  const int id = gp.pdgId();

	  // Only leptons
	  if(abs(id) == 11 || abs(id) == 13 || abs(id) == 15 ||
	     abs(id) == 16 || abs(id) == 14 || abs(id) == 12){

	    const reco::Candidate *it     = 0;
	    const reco::Candidate *mom    = 0;
	    const reco::Candidate *itmom  = 0;
	    const reco::Candidate *mommom = 0;

	    int momid = id;
	    it = (&gp);
	    // This loop searches the particle's mother
	    while(momid == id){
	      if(it != 0){
		mom = it->mother();
		if(mom != 0) momid = mom->pdgId();
		else momid = 0;
		if(momid == id) it = mom;
	      } // if(it != 0)
	      else momid = 0;
	    } // while(momid == id)

	    int mommomid = momid;

	    if(mom != 0){
	      itmom = mom;
	      // This loop searches the mother's mom of the particle
	      while (mommomid == momid){
		if(itmom !=0){
		  mommom = itmom->mother();
		  if(mommom != 0) mommomid = mommom->pdgId();
		  else mommomid = 0;
		  if(mommomid == momid) itmom = mommom->mother();
		}
		else mommomid = 0;
	      } // if(mom != 0)
	    } // while(mommomid == momid)

	    if (abs(momid) == 24 && abs(mommomid) == 6){
	      if (abs(id) == 13 || abs(id) == 11) Genlepton.SetPtEtaPhiE(gp.pt(), gp.eta(), gp.phi(), gp.energy());  

	      if (abs(id) == 16 || abs(id) == 14 || abs(id) == 12) Gennu.SetPtEtaPhiE(gp.pt(), gp.eta(), gp.phi(), gp.energy());  
	      
	      if (abs(id) == 15){ // Taus 
		for(unsigned int h = 0; h <  gp.numberOfDaughters(); h++) {
		  const reco::Candidate *gd = gp.daughter(h);
		  const int taudauid = gd->pdgId();
		  if (abs(taudauid) == 13 || abs(taudauid) == 11) Genlepton.SetPtEtaPhiE(gd->pt(), gd->eta(), gd->phi(), gd->energy());  
		} // for(taus' daughters)
	      } // if(taus)
	        
	    } // if(t->W)
	    
	  }// if (mu || e || tau)
	} // for(genParticles)

	// Full lepton and neutrino information
	b_GenNu_pT  = Gennu.Pt();
	b_GenNu_eta = Gennu.Eta();
	b_GenNu_phi = Gennu.Phi();
	b_GenNu_E   = Gennu.E();
	
	b_GenLepton_pT  = Genlepton.Pt();
	b_GenLepton_eta = Genlepton.Eta();
	b_GenLepton_phi = Genlepton.Phi();
	b_GenLepton_E   = Genlepton.E();

	edm::Handle<reco::GenJetCollection> genJets;
	iEvent.getByToken(genJetToken_, genJets);
	
	edm::Handle<vector<vector<int>>> JetMom;
	iEvent.getByToken(JetMotherToken_, JetMom);	
	
	for (unsigned int j = 0; j < genJets->size(); j++){
	  const cat::GenJet & gjet = (*genJets)[j];
	  if(std::abs(gjet.eta())< 2.5 &&
	     gjet.pt()> 20){
	    
            b_GenJet_pT ->push_back(gjet.pt());
            b_GenJet_eta->push_back(gjet.eta());
            b_GenJet_phi->push_back(gjet.phi());
            b_GenJet_E  ->push_back(gjet.energy());

            b_GenJet_pT->push_back(gjet.pt());

	    std::vector<int> moms = (*JetMom)[j];	    
	    int WIDJet = 0, TopIDJet = 0; 
            for(unsigned int nmom = 0; nmom < moms.size() ; nmom++){
              int momID = moms.at(nmom);
              if(std::abs(momID) == 24) WIDJet   = momID;
              if(std::abs(momID) == 6 ) TopIDJet = momID;
            }

	    if (WIDJet != 0) b_GenJet_mom->push_back(WIDJet);
	    else if (TopIDJet != 0)  b_GenJet_mom->push_back(TopIDJet);
	    else b_GenJet_mom->push_back(0);
 
	    b_GenJet_GenConeMom->push_back(0);
	  }// if(Good genJet)
	}// for(genJet)	

	// GenTop Top, W and Add GEN-Jets 
	int GenConeMomID[6] = {6, 6, 24, 24, 0, 0};
	for (int ijGT=0; ijGT<6; ijGT++){
	  int GenConeGenJetIndex = -999;
	  float mGenConegJet_E = (*b_GenCone_gJet_E)[ijGT];
	  auto itGenConeGenJet = std::find((*b_GenJet_E).begin(), (*b_GenJet_E).end(), mGenConegJet_E);
	  if(itGenConeGenJet != (*b_GenJet_E).end()) GenConeGenJetIndex = std::distance((*b_GenJet_E).begin(), itGenConeGenJet);
	  b_GenCone_gJetIndex->push_back(GenConeGenJetIndex);

	  if(GenConeGenJetIndex != -999) (*b_GenJet_GenConeMom)[GenConeGenJetIndex] = GenConeMomID[ijGT];

	}// for(ijGT)

	if(b_GenJet_pT->size() >= SkimNJets_) gentree->Fill();

      } // if(nGenLep == 1)
    }// if(GENTTbarMCTree_)

  } // if(TTbarMC>0)


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Primary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  edm::Handle<int> pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );

  b_nGoodPV = *pvHandle;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Missing E_T
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  Handle<cat::METCollection> MET;
  iEvent.getByToken(metToken_, MET);

  // MET-PF
  b_MET     = MET->begin()->pt();
  b_MET_phi = MET->begin()->phi();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  cat::ElectronCollection selectedElectrons;
  cat::ElectronCollection vetoElectrons;

  Handle<cat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  for (unsigned int i = 0; i < electrons->size() ; i++) {
    const cat::Electron & electron = electrons->at(i);

    if( IsSelectElectron( electron ) ) selectedElectrons.push_back( electron );
    else if( IsVetoElectron( electron ) ) vetoElectrons.push_back( electron ); // does not Include selected electrons

  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  cat::MuonCollection selectedMuons;
  cat::MuonCollection vetoMuons;

  Handle<cat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (unsigned int i = 0; i < muons->size() ; i++) {
    const cat::Muon & muon = muons->at(i);

    if( IsSelectMuon( muon) ) selectedMuons.push_back( muon);
    else if( IsVetoMuon( muon) ) vetoMuons.push_back( muon); // does not Include selected muons

  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  TLorentzVector lepton;
  std::vector<float> *lepton_SF;
  lepton_SF = new std::vector<float>;

  float lepton_LES = 0.0;
  int ch_tag  = 999;

  if(selectedMuons.size()     == 1 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 0 &&
     vetoElectrons.size()     == 0){
    lepton.SetPtEtaPhiE(selectedMuons[0].pt(), selectedMuons[0].eta(), selectedMuons[0].phi(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets

      if(isMC_) {
	// Lepton SF (ID/ISO)
	lepton_SF->push_back( SF_muon_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()) ) );       // [0]-> SF
	lepton_SF->push_back( SF_muon_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()),  1.0 ) ); // [1]-> SF+Error
	lepton_SF->push_back( SF_muon_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()), -1.0 ) ); // [2]-> SF-Error
	//LES
	lepton_LES = selectedMuons[0].shiftedEn();
      }
  }

  if(selectedMuons.size()     == 0 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 1 &&
     vetoElectrons.size()     == 0){
    ch_tag = 1; //electron + jets
    lepton.SetPtEtaPhiE(selectedElectrons[0].pt(), selectedElectrons[0].eta(), selectedElectrons[0].phi(), selectedElectrons[0].energy());

     if(isMC_) {
       // Lepton SF (ID/ISO)
       lepton_SF->push_back( SF_elec_( selectedElectrons[0].pt(), std::abs(selectedElectrons[0].scEta()) ) );       // [0]-> SF
       lepton_SF->push_back( SF_elec_( selectedElectrons[0].pt(), std::abs(selectedElectrons[0].scEta()),  1.0 ) ); // [1]-> SF+Error
       lepton_SF->push_back( SF_elec_( selectedElectrons[0].pt(), std::abs(selectedElectrons[0].scEta()), -1.0 ) ); // [2]-> SF-Error
       //LES
       lepton_LES = selectedElectrons[0].shiftedEn();
     }
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // HLTrigger
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  std::vector<string> triggerNameEl_, triggerNameMu_;
  if(isMC_){
    triggerNameEl_ = triggerNameMCEl_;
    triggerNameMu_ = triggerNameMCMu_;
  } 
  else{
    triggerNameEl_ = triggerNameDataEl_;
    triggerNameMu_ = triggerNameDataMu_;  
  }  
  
  bool EvTrigger = false; 
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  if( !iEvent.getByToken(triggerBits_, triggerBits) ){
    iEvent.getByToken(triggerBits2_, triggerBits);
  } 
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  bool IsTriggerMu = false;
  bool IsTriggerEl = false;
  
  // DEBUG: Print all triggers
  // trigHelper.listFiredTriggers();

  for(std::vector<string>::iterator TrMu_it = triggerNameMu_.begin(); TrMu_it != triggerNameMu_.end(); TrMu_it++){
    // No trigger
    IsTriggerMu = trigHelper.triggerFired(*TrMu_it);
    if(*TrMu_it == "notrigger") IsTriggerMu = true; 
    if (IsTriggerMu) break;
  }
  
  for(std::vector<string>::iterator TrEl_it = triggerNameEl_.begin(); TrEl_it != triggerNameEl_.end(); TrEl_it++){
    IsTriggerEl = trigHelper.triggerFired(*TrEl_it);
    if(*TrEl_it == "notrigger") IsTriggerEl = true;
    if (IsTriggerEl) break;
  }

  if ( (ch_tag == 0 && IsTriggerMu) ||
       (ch_tag == 1 && IsTriggerEl) 
       ) {
    EvTrigger = true;
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Fill Tree with events that have ONLY one lepton
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Check Gen Level for ttbar sample
  if (TTbarMC_ >0){
    if(TTbarMC_ == 1){ // Signal ttbar event
      if(nGenLep != 1) ch_tag = 999;
    }
    if(TTbarMC_ == 2){ // Background ttbar event
      if(nGenLep == 1) ch_tag = 999;
    }

    // ttbar Category
    if(!IsCat) ch_tag = 999;

  } // if(TTbarMC_ >0)

  if (ch_tag<2 && EvTrigger){ // Single lepton event

    b_Channel  = ch_tag;

    b_Lepton_pT  = lepton.Pt();
    b_Lepton_eta = lepton.Eta();
    b_Lepton_phi = lepton.Phi();
    b_Lepton_E   = lepton.E();

    b_Lepton_SF  = lepton_SF;
    b_Lepton_LES = lepton_LES;

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Jets
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    Handle<cat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    // Arrange Jets by CSV discriminator
    std::vector<int>  JetIndex;
    for (unsigned int i = 0; i < jets->size() ; i++) JetIndex.push_back(i);
    for(unsigned int ijet=0; ijet < JetIndex.size(); ijet++){
      for(unsigned int jjet=ijet+1; jjet < JetIndex.size(); jjet++){
        const cat::Jet & jet_i = jets->at(JetIndex[ijet]);
        const cat::Jet & jet_j = jets->at(JetIndex[jjet]);

        float iJetCSV = jet_i.bDiscriminator(BTAG_CSVv2);
        float jJetCSV = jet_j.bDiscriminator(BTAG_CSVv2);

        if(jJetCSV > iJetCSV){
          float tempIndex = JetIndex[ijet];
          JetIndex[ijet] = JetIndex[jjet];
          JetIndex[jjet] = tempIndex;
        }//if(jJetCSV > iJetCSV)
      }
    }

    int N_GoodJets = 0;

    // Initialize SF_btag
    float Jet_SF_CSV[4][19];
    for (unsigned int iu=0; iu<19; iu++){
      for (unsigned int ipTj=0; ipTj<4; ipTj++) Jet_SF_CSV[ipTj][iu] = 1.0;
    }
    // Run again over all Jets (CSV order)
    for (unsigned int i = 0; i < JetIndex.size() ; i++) {

      const cat::Jet & jet = jets->at(JetIndex[i]);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection (pT>15GeV to take into account SYST Variations)
      if(std::abs(jet.eta()) < 2.4 && jet.pt() > 15. && jet.LooseId()) goodJet = true;
      // Jet Cleaning
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if(dr_LepJet > 0.4) cleanJet = true;

      if(goodJet && cleanJet){
        // Basic variables
        N_GoodJets ++;
        b_Jet_pT ->push_back(jet.pt());
        b_Jet_eta->push_back(jet.eta());
        b_Jet_phi->push_back(jet.phi());
        b_Jet_E  ->push_back(jet.energy());
        b_Jet_Index ->push_back(JetIndex[i]);

        // Parton Flavour
        b_Jet_partonFlavour->push_back(jet.partonFlavour());
        b_Jet_hadronFlavour->push_back(jet.hadronFlavour());

        // b-tag discriminant
        float jet_btagDis_CSV = jet.bDiscriminator(BTAG_CSVv2);
        b_Jet_CSV ->push_back(jet_btagDis_CSV);
        // c-tag discriminant
        float jet_btagDis_CvsL = jet.bDiscriminator("pfCombinedCvsLJetTags");
        b_Jet_CvsL ->push_back(jet_btagDis_CvsL);
        float jet_btagDis_CvsB = jet.bDiscriminator("pfCombinedCvsBJetTags");
	b_Jet_CvsB ->push_back(jet_btagDis_CvsB);

        if(isMC_) {
          // JES
          b_Jet_JES_Up  ->push_back(jet.shiftedEnUp());
          b_Jet_JES_Down->push_back(jet.shiftedEnDown());

          // JER
          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
          b_Jet_JER_Up   ->push_back(jet.smearedResUp());
          b_Jet_JER_Nom  ->push_back(jet.smearedRes());
          b_Jet_JER_Down ->push_back(jet.smearedResDown());

          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
          // Saving the central SF and the 18 syst. unc. for pT_Jets > 25 GeV
          if(jet.pt() > 25.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[0][iu] *= SF_CSV_.getSF(jet, iu);
          // Saving the central SF and the 18 syst. unc. for pT_Jets > 30 GeV
          if(jet.pt() > 30.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[1][iu] *= SF_CSV_.getSF(jet, iu);
          // Saving the central SF and the 18 syst. unc. for pT_Jets > 35 GeV
          if(jet.pt() > 35.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[2][iu] *= SF_CSV_.getSF(jet, iu);
          // Saving the central SF and the 18 syst. unc. for pT_Jets > 40 GeV
          if(jet.pt() > 40.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[3][iu] *= SF_CSV_.getSF(jet, iu);

	  // GEN-Jets matched with RECO-Jets
	  // Reference to gen object
	  if(TTbarMC_== 1){
	    int MatchedGenJetIndex = -999;
	    int rJetGenConeMom = -999;
	    if(jet.genJet()){ 
	      float mGenJet_E = jet.genJet()->energy();
	      auto itGenJet = std::find((*b_GenJet_E).begin(), (*b_GenJet_E).end(), mGenJet_E);
	      if(itGenJet != (*b_GenJet_E).end()) MatchedGenJetIndex = std::distance((*b_GenJet_E).begin(), itGenJet);

	      if(MatchedGenJetIndex != -999){
		if      (MatchedGenJetIndex == (*b_GenCone_gJetIndex)[0] ||
		         MatchedGenJetIndex == (*b_GenCone_gJetIndex)[1]) rJetGenConeMom = 6;
		else if (MatchedGenJetIndex == (*b_GenCone_gJetIndex)[2] ||
			 MatchedGenJetIndex == (*b_GenCone_gJetIndex)[3]) rJetGenConeMom = 24;
		else if (MatchedGenJetIndex == (*b_GenCone_gJetIndex)[4] ||
			 MatchedGenJetIndex == (*b_GenCone_gJetIndex)[5]) rJetGenConeMom = 0;
	      } // if(MatchedGenJetIndex != -999) 
	    } // if(jet.genJet())

	    b_Jet_MatchedGenJetIndex->push_back(MatchedGenJetIndex);

	    b_Jet_GenConeMom->push_back(rJetGenConeMom);
	    
	  } // if(TTbarMC_== 1)
        } // if(isMC_)	
      }// if(GoodJets)
    }// for(AllJets)
    
    b_Jet_Number = N_GoodJets;

    for (unsigned int iu=0; iu<19; iu++){
      b_Jet_SF_CSV_25->push_back(1.0);
      b_Jet_SF_CSV_30->push_back(1.0);
      b_Jet_SF_CSV_35->push_back(1.0);
      b_Jet_SF_CSV_40->push_back(1.0);
    }
    (*b_Jet_SF_CSV_25)[0] = Jet_SF_CSV[0][0]; //Central
    (*b_Jet_SF_CSV_30)[0] = Jet_SF_CSV[1][0]; //Central
    (*b_Jet_SF_CSV_35)[0] = Jet_SF_CSV[2][0]; //Central
    (*b_Jet_SF_CSV_40)[0] = Jet_SF_CSV[3][0]; //Central
    // To save only the error
    for (unsigned int iu=1; iu<19; iu++){
      (*b_Jet_SF_CSV_25)[iu] = std::abs(Jet_SF_CSV[0][iu] - Jet_SF_CSV[0][0]) ; // Syst. Unc.
      (*b_Jet_SF_CSV_30)[iu] = std::abs(Jet_SF_CSV[1][iu] - Jet_SF_CSV[1][0]) ; // Syst. Unc.
      (*b_Jet_SF_CSV_35)[iu] = std::abs(Jet_SF_CSV[2][iu] - Jet_SF_CSV[2][0]) ; // Syst. Unc.
      (*b_Jet_SF_CSV_40)[iu] = std::abs(Jet_SF_CSV[3][iu] - Jet_SF_CSV[3][0]) ; // Syst. Unc.
    }
    
    //---------------------------------------------------------------------------
    // Kinematic Reconstruction
    //---------------------------------------------------------------------------
    TLorentzVector Kinnu, Kinblrefit, Kinbjrefit, Kinj1refit, Kinj2refit;
    Kinnu.SetPtEtaPhiE(0,0,0,0);
    Kinblrefit.SetPtEtaPhiE(0,0,0,0);
    Kinbjrefit.SetPtEtaPhiE(0,0,0,0);
    Kinj1refit.SetPtEtaPhiE(0,0,0,0);
    Kinj2refit.SetPtEtaPhiE(0,0,0,0);
    
    std::vector<int> KinBestIndices;
    KinBestIndices.push_back(-999);
    KinBestIndices.push_back(-999);
    KinBestIndices.push_back(-999);
    KinBestIndices.push_back(-999);
    float bestchi2 = 0;
    
    if(N_GoodJets > 3){
      
      TLorentzVector KinLep;
      KinLep = lepton;
      
      TLorentzVector KinMET;
      KinMET.SetPtEtaPhiE(b_MET, 0.0, b_MET_phi, b_MET);
      
      std::vector<ComJet> KinJets;
      for (unsigned int kj=0; kj<b_Jet_pT->size(); kj++){
	ComJet kjet;
	kjet.SetPtEtaPhiE((*b_Jet_pT)[kj],(*b_Jet_eta)[kj],(*b_Jet_phi)[kj],(*b_Jet_E)[kj]);
	kjet.CSV = (*b_Jet_CSV)[kj];
	
	KinJets.push_back(kjet);
      }
      
      FindHadronicTop(KinLep, KinJets, KinMET, KFUsebtag_, CSVPosConKF_, KinBestIndices, bestchi2, Kinnu, Kinblrefit, Kinbjrefit, Kinj1refit, Kinj2refit);
      
      // for (unsigned int iin =0; iin<KinBestIndices.size(); iin++) std::cout << KinBestIndices.at(iin) << std::endl;
      // std::cout << "Best Chi2 = " << bestchi2 << std::endl;
      // std::cout << "Lep pT (NEW) = " << KinLep.Pt() << std::endl;
      
      // std::cout << "\nJet[0] pT = " << KinJets[0].Pt() << std::endl;
      // std::cout << "Jet[1] pT = " << KinJets[1].Pt() << std::endl;
      // std::cout << "Jet[2] pT = " << KinJets[2].Pt() << std::endl;
      // std::cout << "Jet[3] pT = " << KinJets[3].Pt() << std::endl;
      // std::cout << "Jet[4] pT = " << KinJets[4].Pt() << std::endl;
      // std::cout << "Jet[5] pT = " << KinJets[5].Pt() << std::endl;
      
      // std::cout << "\nMET = " << KinMET.Et() << std::endl;
      // std::cout << "nu pT = " << Kinnu.Pt() << std::endl;
      // std::cout << "nu ET = " << Kinnu.Et() << std::endl;
      // std::cout << "nu E = " << Kinnu.E() << std::endl;
    
    } //if(N_GoodJets > 3)   
    
    // LorentzVector for Jets. Same order as KinBestIndices
    b_KinJet_pT->push_back(Kinbjrefit.Pt());
    b_KinJet_pT->push_back(Kinj1refit.Pt());
    b_KinJet_pT->push_back(Kinj2refit.Pt()); 
    b_KinJet_pT->push_back(Kinblrefit.Pt());
    
    b_KinJet_eta->push_back(Kinbjrefit.Eta());
    b_KinJet_eta->push_back(Kinj1refit.Eta());
    b_KinJet_eta->push_back(Kinj2refit.Eta()); 
    b_KinJet_eta->push_back(Kinblrefit.Eta());
    
    b_KinJet_phi->push_back(Kinbjrefit.Phi());
    b_KinJet_phi->push_back(Kinj1refit.Phi());
    b_KinJet_phi->push_back(Kinj2refit.Phi()); 
    b_KinJet_phi->push_back(Kinblrefit.Phi());
    
    b_KinJet_E->push_back(Kinbjrefit.E());
    b_KinJet_E->push_back(Kinj1refit.E());
    b_KinJet_E->push_back(Kinj2refit.E()); 
    b_KinJet_E->push_back(Kinblrefit.E());
    
    b_KinNu_pT  = Kinnu.Pt();
    b_KinNu_eta = Kinnu.Eta();
    b_KinNu_phi = Kinnu.Phi();
    b_KinNu_E   = Kinnu.E();
    
    for(unsigned int iki=0; iki<KinBestIndices.size(); iki++) b_KinJet_Index->push_back(KinBestIndices.at(iki));
    
    b_Kin_Chi2 = bestchi2;
        
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Fill Tree with event at 1 lepton cut level
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    
    if(b_Jet_pT->size() >= SkimNJets_) tree->Fill();
        
  }// if(ch_tag)
  
}

//------------- Good Muon Selection -----------------------
bool ttbbLepJetsAnalyzer::IsSelectMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Tight selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.isTightMuon());

  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 20);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.1); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.15 );

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodMuon;
}
//------------- Loose Muon Selection -----------------------
bool ttbbLepJetsAnalyzer::IsVetoMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Loose selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.isLooseMuon());

  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 15);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.25 );

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool ttbbLepJetsAnalyzer::IsSelectElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.1);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // Electron cut based selection
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Summer16-80X-V1-medium") > 0.0;

  // Electron MVA selection (Tight: WP80)
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_for_7_4_12_Spring15_MVA
  // if (std::abs(i_electron_candidate.eta()) < 0.8)        GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.988153; // Barrel
  // else if (std::abs(i_electron_candidate.eta()) > 0.8 &&
  // 	   std::abs(i_electron_candidate.eta()) < 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.967910; // Barrel
  // else if (std::abs(i_electron_candidate.eta()) > 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.841729; // EndCaps

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;

  if ( std::abs(i_electron_candidate.scEta()) <= 1.479)   GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.0766 );
  else GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.0678 );


  // Effective Area Parametrization can be found in:
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}
//------------- Loose Electron Selection -----------------------
bool ttbbLepJetsAnalyzer::IsVetoElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Summer16-80X-V1-veto") > 0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  if ( std::abs(i_electron_candidate.scEta()) <= 1.479) GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.126 );
  else GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.144 );
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbbLepJetsAnalyzer);

