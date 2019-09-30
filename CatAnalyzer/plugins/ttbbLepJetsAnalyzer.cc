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
// Original Author:  Javier Brochero Cifuentes
// Author : Seohyeon An
//         Created:  Tue Feb  3 09:52:55 CET 2015

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
  float weight_valid   (float weight);

  bool isMC_ ;
  bool doLooseLepton_;

  int TTbarMC_; // 0->No ttbar, 1->ttbar Signal, 2->ttbar Background
  int TTbarCatMC_;
    
  // Event Weights
  edm::EDGetTokenT<float>                        genWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           pdfWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           scaleUpWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           scaleDownWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           psWeightToken_;
  edm::EDGetTokenT<float>                        puWeightToken_;
  edm::EDGetTokenT<float>                        puUpWeightToken_;
  edm::EDGetTokenT<float>                        puDownWeightToken_;
  // Object Collections
  edm::EDGetTokenT<int>                          genttbarHiggsCatToken_;
  edm::EDGetTokenT<cat::GenTopCollection>        genttbarCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
  edm::EDGetTokenT<cat::JetCollection>           jetToken_;
  edm::EDGetTokenT<cat::METCollection>           metToken_;
  edm::EDGetTokenT<int>                          pvToken_;
  edm::EDGetTokenT<int>                          nTrueVertToken_;
  edm::EDGetTokenT<int>                          recoFiltersToken_;
  // Trigger
  edm::EDGetTokenT<int>                          trigMuFiltersToken_;
  edm::EDGetTokenT<int>                          trigElFiltersToken_;
  edm::EDGetTokenT<int>                          trigElHTFiltersToken_;

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
  std::vector<float> *b_PDFWeight;
  std::vector<float> *b_PSWeight;
  // PU/Vertices
  std::vector<float> *b_PUWeight;
  int b_nGoodPV, b_nTruePV;
  // Channel and Categorization
  int b_GenChannel;
  int b_Channel;
  std::vector<int>   *b_GenConeCatID;
  std::vector<float> *b_GenCone_gJet_pt;
  std::vector<float> *b_GenCone_gJet_eta;
  std::vector<float> *b_GenCone_gJet_phi;
  std::vector<float> *b_GenCone_gJet_e;
  std::vector<int>   *b_GenCone_gJetFlavW;
  int b_GenCone_NgJetsW;
  int b_GenHiggsCatID;
  // MET
  float b_MET, b_MET_phi;
  // GEN Leptons
  float b_GenLepton_pt;
  float b_GenLepton_eta;
  float b_GenLepton_phi;
  float b_GenLepton_e;
  // GEN Neutrino
  float b_GenNu_pt;
  float b_GenNu_eta;
  float b_GenNu_phi;
  float b_GenNu_e;

  // Leptons
  float b_Lepton_pt;
  float b_Lepton_eta;
  float b_Lepton_phi;
  float b_Lepton_e;
  float b_Lepton_relIso;
  bool b_Lepton_isIso;
  std::vector<float> *b_Lepton_SF;

  // additional b jets 
  float b_addbjet1_pt; 
  float b_addbjet1_eta; 
  float b_addbjet1_phi; 
  float b_addbjet1_e; 

  float b_addbjet2_pt;
  float b_addbjet2_eta;
  float b_addbjet2_phi;
  float b_addbjet2_e; 

  // minimum deltaR b jets
  float b_mindRjet1_pt;
  float b_mindRjet1_eta;
  float b_mindRjet1_phi;
  float b_mindRjet1_e;

  float b_mindRjet2_pt;
  float b_mindRjet2_eta;
  float b_mindRjet2_phi;
  float b_mindRjet2_e;

  float b_mindR;

  // Jets
  int b_Jet_Number;
  std::vector<float> *b_Jet_pt;
  std::vector<float> *b_Jet_eta;
  std::vector<float> *b_Jet_phi;
  std::vector<float> *b_Jet_e;
  std::vector<int>   *b_Jet_Index;
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
  std::vector<float> *b_Jet_deepCSV;
  std::vector<float> *b_Jet_deepJet;
  std::vector<float> *b_Jet_SF_deepCSV_30;
  std::vector<float> *b_Jet_SF_deepJet_30;
  // c-Jet discriminant
  std::vector<float> *b_Jet_deepCvsL,    *b_Jet_deepCvsB;
  std::vector<float> *b_Jet_deepJetCvsL, *b_Jet_deepJetCvsB;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Histograms
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Histograms: Number of Events and Weights
  TH1D *EventInfo, *ScaleWeights, *PDFWeights, *PSWeights;
  // Scale factor evaluators
  BTagWeightEvaluator SF_deepCSV_;
  ScaleFactorEvaluator SF_muonId_, SF_muonIso_, SF_muonTrg_,
                       SF_elecId_, SF_elecReco_, SF_elecZvtx_, SF_elecTrg_;
 
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
  doLooseLepton_(iConfig.getUntrackedParameter<bool>("doLooseLepton", false)),
  TTbarMC_    (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0)),
  TTbarCatMC_ (iConfig.getUntrackedParameter<int>("TTbarCatLabel", 0))
{
  const auto elecIdSFSet = iConfig.getParameter<edm::ParameterSet>("elecIdSF");
  SF_elecId_.set(elecIdSFSet.getParameter<std::vector<double>>("pt_bins" ),
	         elecIdSFSet.getParameter<std::vector<double>>("eta_bins"),
	         elecIdSFSet.getParameter<std::vector<double>>("values"  ),
	         elecIdSFSet.getParameter<std::vector<double>>("errors"  ));
  const auto elecRecoSFSet = iConfig.getParameter<edm::ParameterSet>("elecRecoSF");
  SF_elecReco_.set(elecRecoSFSet.getParameter<std::vector<double>>("pt_bins" ),
	           elecRecoSFSet.getParameter<std::vector<double>>("eta_bins"),
	           elecRecoSFSet.getParameter<std::vector<double>>("values"  ),
	           elecRecoSFSet.getParameter<std::vector<double>>("errors"  ));
  const auto elecZvtxSet = iConfig.getParameter<edm::ParameterSet>("elecZvtxSF");
  SF_elecZvtx_.set(elecZvtxSet.getParameter<std::vector<double>>("pt_bins" ),
	           elecZvtxSet.getParameter<std::vector<double>>("eta_bins"),
	           elecZvtxSet.getParameter<std::vector<double>>("values"  ),
	           elecZvtxSet.getParameter<std::vector<double>>("errors"  ));
  const auto elecTrgSFSet = iConfig.getParameter<edm::ParameterSet>("elecTrgSF");
  SF_elecTrg_.set(elecTrgSFSet.getParameter<std::vector<double>>("pt_bins" ),
	          elecTrgSFSet.getParameter<std::vector<double>>("eta_bins"),
	          elecTrgSFSet.getParameter<std::vector<double>>("values"  ),
	          elecTrgSFSet.getParameter<std::vector<double>>("errors"  ));

  const auto muonIdSFSet = iConfig.getParameter<edm::ParameterSet>("muonIdSF");
  SF_muonId_.set(muonIdSFSet.getParameter<std::vector<double>>("pt_bins"    ),
	         muonIdSFSet.getParameter<std::vector<double>>("abseta_bins"),
	         muonIdSFSet.getParameter<std::vector<double>>("values"     ),
	         muonIdSFSet.getParameter<std::vector<double>>("errors"     ));
  const auto muonIsoSFSet = iConfig.getParameter<edm::ParameterSet>("muonIsoSF");
  SF_muonIso_.set(muonIsoSFSet.getParameter<std::vector<double>>("pt_bins"    ),
	       muonIsoSFSet.getParameter<std::vector<double>>("abseta_bins"),
	       muonIsoSFSet.getParameter<std::vector<double>>("values"     ),
	       muonIsoSFSet.getParameter<std::vector<double>>("errors"     ));
  const auto muonTrgSFSet = iConfig.getParameter<edm::ParameterSet>("muonTrgSF");
  SF_muonTrg_.set(muonTrgSFSet.getParameter<std::vector<double>>("pt_bins"    ),
	          muonTrgSFSet.getParameter<std::vector<double>>("abseta_bins"),
	          muonTrgSFSet.getParameter<std::vector<double>>("values"     ),
	          muonTrgSFSet.getParameter<std::vector<double>>("errors"     ));
 
  SF_deepCSV_.initCSVWeight(false, "deepcsv");
  
  // Weights
  auto genWeightLabel = iConfig.getParameter<edm::InputTag>("genWeightLabel");
  // aMC@NLO
  genWeightToken_       = consumes<float>             (edm::InputTag(genWeightLabel.label()));
  // PDF
  pdfWeightToken_       = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "pdf"));
  // Scale
  scaleUpWeightToken_   = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaleup"));
  scaleDownWeightToken_ = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaledown"));
  // PS
  psWeightToken_        = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "ps"));
  // PileUp
  auto puWeightLabel    = iConfig.getParameter<edm::InputTag>("puWeightLabel");
  puWeightToken_        = consumes<float>             (edm::InputTag(puWeightLabel.label()));
  puUpWeightToken_      = consumes<float>             (edm::InputTag(puWeightLabel.label(),"up"));
  puDownWeightToken_    = consumes<float>             (edm::InputTag(puWeightLabel.label(),"dn"));
  // CSV Weights
  auto deepcsvWeightLabel = iConfig.getParameter<edm::InputTag>("deepcsvWeightLabel");
  // GEN and ttbar Categorization
  genttbarCatToken_      = consumes<cat::GenTopCollection>        (iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));
  genttbarHiggsCatToken_ = consumes<int>                          (iConfig.getParameter<edm::InputTag>("genHiggsCatLabel"));
  // Object Collections
  muonToken_         = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_     = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_          = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_          = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_           = consumes<int>                          (iConfig.getParameter<edm::InputTag>("pvLabel"));
  // Trigger 
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  trigMuFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMuFilters"));
  trigElFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigElFilters"));
  trigElHTFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigElHTFilters"));
  // PU = 0 prescription
  nTrueVertToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nTrueVertLabel"));

  b_PUWeight      = new std::vector<float>;
  b_PDFWeight     = new std::vector<float>;  
  b_ScaleWeight   = new std::vector<float>;  
  b_Lepton_SF     = new std::vector<float>;
  b_PSWeight      = new std::vector<float>;

  b_GenConeCatID      = new std::vector<int>;
  b_GenCone_gJet_pt   = new std::vector<float>;
  b_GenCone_gJet_eta  = new std::vector<float>;
  b_GenCone_gJet_phi  = new std::vector<float>;
  b_GenCone_gJet_e    = new std::vector<float>;
  b_GenCone_gJetFlavW = new std::vector<int>;

  b_Jet_pt            = new std::vector<float>;
  b_Jet_eta           = new std::vector<float>;
  b_Jet_phi           = new std::vector<float>;
  b_Jet_e             = new std::vector<float>;
  b_Jet_Index         = new std::vector<int>;
  b_Jet_deepCSV       = new std::vector<float>;
  b_Jet_deepJet       = new std::vector<float>;
  b_Jet_SF_deepCSV_30 = new std::vector<float>;
  b_Jet_SF_deepJet_30 = new std::vector<float>;
  b_Jet_deepCvsL      = new std::vector<float>;  
  b_Jet_deepCvsB      = new std::vector<float>;  
  b_Jet_deepJetCvsL   = new std::vector<float>;
  b_Jet_deepJetCvsB   = new std::vector<float>;
  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;
  b_Jet_JES_Up        = new std::vector<float>;
  b_Jet_JES_Down      = new std::vector<float>;
  b_Jet_JER_Up        = new std::vector<float>;
  b_Jet_JER_Nom       = new std::vector<float>;
  b_Jet_JER_Down      = new std::vector<float>;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "event/I");
  tree->Branch("run",        &b_Run,         "run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "luminumber/I");
  tree->Branch("genweight",  &b_GenWeight,   "genweight/F");
  tree->Branch("TruePV",     &b_nTruePV,     "TruePV/I");
  tree->Branch("GoodPV",     &b_nGoodPV,     "GoodPV/I");
  tree->Branch("channel",    &b_Channel,     "channel/I");

  tree->Branch("PUWeight",      "std::vector<float>", &b_PUWeight);
  tree->Branch("pdfweight",     "std::vector<float>", &b_PDFWeight);
  tree->Branch("scaleweight",   "std::vector<float>", &b_ScaleWeight);
  tree->Branch("psweight",      "std::vector<float>", &b_PSWeight);

  tree->Branch("MET",     &b_MET,     "MET/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("lepton_pt",  &b_Lepton_pt,  "lepton_pt/F");
  tree->Branch("lepton_eta", &b_Lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi", &b_Lepton_phi, "lepton_phi/F");
  tree->Branch("lepton_e" ,  &b_Lepton_e,   "lepton_e/F" );

  tree->Branch("lepton_SF",  "std::vector<float>", &b_Lepton_SF );

  tree->Branch("lepton_relIso", &b_Lepton_relIso, "lepton_relIso/F");
  tree->Branch("lepton_isIso",  &b_Lepton_isIso,  "lepton_isIso/O");

  tree->Branch("jet_pt",            "std::vector<float>", &b_Jet_pt);
  tree->Branch("jet_eta",           "std::vector<float>", &b_Jet_eta);
  tree->Branch("jet_phi",           "std::vector<float>", &b_Jet_phi);
  tree->Branch("jet_e" ,            "std::vector<float>", &b_Jet_e);
  tree->Branch("jet_index" ,        "std::vector<int>",   &b_Jet_Index);
  tree->Branch("jet_deepCSV" ,      "std::vector<float>", &b_Jet_deepCSV);
  tree->Branch("jet_deepJet" ,      "std::vector<float>", &b_Jet_deepJet);
  tree->Branch("jet_SF_deepCSV_30", "std::vector<float>", &b_Jet_SF_deepCSV_30);
  tree->Branch("jet_SF_deepJet_30", "std::vector<float>", &b_Jet_SF_deepJet_30);
  tree->Branch("jet_deepCvsL",      "std::vector<float>", &b_Jet_deepCvsL);
  tree->Branch("jet_deepCvsB",      "std::vector<float>", &b_Jet_deepCvsB);
  tree->Branch("jet_deepJetCvsL",   "std::vector<float>", &b_Jet_deepJetCvsL);
  tree->Branch("jet_deepJetCvsB",   "std::vector<float>", &b_Jet_deepJetCvsB);

  tree->Branch("jet_number" , &b_Jet_Number, "jet_number/I" );

  tree->Branch("jet_partonFlavour", "std::vector<int>",   &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>",   &b_Jet_hadronFlavour);
  tree->Branch("jet_JES_Up",        "std::vector<float>", &b_Jet_JES_Up);
  tree->Branch("jet_JES_Down",      "std::vector<float>", &b_Jet_JES_Down);
  tree->Branch("jet_JER_Up",        "std::vector<float>", &b_Jet_JER_Up);
  tree->Branch("jet_JER_Nom",       "std::vector<float>", &b_Jet_JER_Nom);
  tree->Branch("jet_JER_Down",      "std::vector<float>", &b_Jet_JER_Down);

  // GEN Variables (only ttbarSignal)
  if(TTbarMC_ == 1){
    tree->Branch("pdfweight",   "std::vector<float>", &b_PDFWeight);
    tree->Branch("scaleweight", "std::vector<float>", &b_ScaleWeight);

    tree->Branch("genconecatid" ,      "std::vector<int>",   &b_GenConeCatID);
    tree->Branch("gencone_gjet_pt" ,   "std::vector<float>", &b_GenCone_gJet_pt);
    tree->Branch("gencone_gjet_eta" ,  "std::vector<float>", &b_GenCone_gJet_eta);
    tree->Branch("gencone_gjet_phi" ,  "std::vector<float>", &b_GenCone_gJet_phi);
    tree->Branch("gencone_gjet_e" ,    "std::vector<float>", &b_GenCone_gJet_e);
    tree->Branch("gencone_gJetFlavW" , "std::vector<int>",   &b_GenCone_gJetFlavW);

    tree->Branch("gencone_NgjetsW", &b_GenCone_NgJetsW, "gencone_NgjetsW/I");
    tree->Branch("draddjets",       &b_DRAddJets,       "draddjets/F");
    tree->Branch("genhiggscatid",   &b_GenHiggsCatID,   "genhiggscatid/I");

    tree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    tree->Branch("genlepton_pt",  &b_GenLepton_pt,  "genlepton_pt/F");
    tree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    tree->Branch("genlepton_phi", &b_GenLepton_phi, "genlepton_phi/F");
    tree->Branch("genlepton_e",   &b_GenLepton_e,   "genlepton_e/F");
    tree->Branch("gennu_pt",      &b_GenNu_pt,      "gennu_pt/F");
    tree->Branch("gennu_eta",     &b_GenNu_eta,     "gennu_eta/F");
    tree->Branch("gennu_phi",     &b_GenNu_phi,     "gennu_phi/F");
    tree->Branch("gennu_e",       &b_GenNu_e,       "gennu_e/F");


    tree->Branch("addbjet1_pt",  &b_addbjet1_pt,  "addbjet1_pt/F"); 
    tree->Branch("addbjet1_eta", &b_addbjet1_eta, "addbjet1_eta/F"); 
    tree->Branch("addbjet1_phi", &b_addbjet1_phi, "addbjet1_phi/F"); 
    tree->Branch("addbjet1_e",   &b_addbjet1_e,   "addbjet1_e/F"); 
    tree->Branch("addbjet2_pt",  &b_addbjet2_pt,  "addbjet2_pt/F");
    tree->Branch("addbjet2_eta", &b_addbjet2_eta, "addbjet2_eta/F");
    tree->Branch("addbjet2_phi", &b_addbjet2_phi, "addbjet2_phi/F");
    tree->Branch("addbjet2_e",   &b_addbjet2_e,   "addbjet2_e/F");  

    tree->Branch("mindRjet1_pt",  &b_mindRjet1_pt,  "mindRjet1_pt/F");
    tree->Branch("mindRjet1_eta", &b_mindRjet1_eta, "mindRjet1_eta/F");
    tree->Branch("mindRjet1_phi", &b_mindRjet1_phi, "mindRjet1_phi/F");
    tree->Branch("mindRjet1_e",   &b_mindRjet1_e,   "mindRjet1_e/F");
    tree->Branch("mindRjet2_pt",  &b_mindRjet2_pt,  "mindRjet2_pt/F");
    tree->Branch("mindRjet2_eta", &b_mindRjet2_eta, "mindRjet2_eta/F");
    tree->Branch("mindRjet2_phi", &b_mindRjet2_phi, "mindRjet2_phi/F");
    tree->Branch("mindRjet2_e",   &b_mindRjet2_e,   "mindRjet2_e/F");
    tree->Branch("mindR",         &b_mindR,         "mimdR/F");

    //GEN TREE
    gentree = fs->make<TTree>("gentree", "TopGENTree");
    gentree->Branch("genweight",     &b_GenWeight,     "genweight/F");
    gentree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    gentree->Branch("genhiggscatid", &b_GenHiggsCatID, "genhiggscatid/I");
    gentree->Branch("draddjets",     &b_DRAddJets,     "draddjets/F");
    gentree->Branch("genlepton_pt",  &b_GenLepton_pt,  "genlepton_pt/F");
    gentree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");

    gentree->Branch("scaleweight",        "std::vector<float>", &b_ScaleWeight );
    gentree->Branch("genconecatid" ,      "std::vector<int>",   &b_GenConeCatID);

    gentree->Branch("addbjet1_pt",  &b_addbjet1_pt,  "addbjet1_pt/F"); 
    gentree->Branch("addbjet1_eta", &b_addbjet1_eta, "addbjet1_eta/F"); 
    gentree->Branch("addbjet1_phi", &b_addbjet1_phi, "addbjet1_phi/F"); 
    gentree->Branch("addbjet1_e",   &b_addbjet1_e,   "addbjet1_e/F"); 
    gentree->Branch("addbjet2_pt",  &b_addbjet2_pt,  "addbjet2_pt/F");
    gentree->Branch("addbjet2_eta", &b_addbjet2_eta, "addbjet2_eta/F");
    gentree->Branch("addbjet2_phi", &b_addbjet2_phi, "addbjet2_phi/F");
    gentree->Branch("addbjet2_e",   &b_addbjet2_e,   "addbjet2_e/F");  

    gentree->Branch("mindRjet1_pt",  &b_mindRjet1_pt,  "mindRjet1_pt/F");
    gentree->Branch("mindRjet1_eta", &b_mindRjet1_eta, "mindRjet1_eta/F");
    gentree->Branch("mindRjet1_phi", &b_mindRjet1_phi, "mindRjet1_phi/F");
    gentree->Branch("mindRjet1_e",   &b_mindRjet1_e,   "mindRjet1_e/F");
    gentree->Branch("mindRjet2_pt",  &b_mindRjet2_pt,  "mindRjet2_pt/F");
    gentree->Branch("mindRjet2_eta", &b_mindRjet2_eta, "mindRjet2_eta/F");
    gentree->Branch("mindRjet2_phi", &b_mindRjet2_phi, "mindRjet2_phi/F");
    gentree->Branch("mindRjet2_e",   &b_mindRjet2_e,   "mindRjet2_e/F");
    gentree->Branch("mindR",         &b_mindR,         "mimdR/F");

  }

  EventInfo = fs->make<TH1D>("EventInfo","Event Information",9,0,9);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
  EventInfo->GetXaxis()->SetBinLabel(3,"ttbb Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(4,"ttbb Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(5,"ttbb Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(6,"ttjj Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(7,"ttjj Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(8,"ttjj Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(9,"Sum of PU Weights");

  ScaleWeights = fs->make<TH1D>("ScaleWeights","Event Weights:",6,0,6);
  ScaleWeights->GetXaxis()->SetBinLabel(1,"muR=Nom  muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(2,"muR=Nom  muF=Down");
  ScaleWeights->GetXaxis()->SetBinLabel(3,"muR=Up   muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(4,"muR=Up   muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(5,"muR=Down muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(6,"muR=Down muF=Down");

  PSWeights = fs->make<TH1D>("PSWeights", "PS Weights",4,0,4);
  PSWeights->GetXaxis()->SetBinLabel(1, "isrDefHi isr:muRfac=0.5");
  PSWeights->GetXaxis()->SetBinLabel(2, "fsrDefHi fsr:muRfac=0.5");
  PSWeights->GetXaxis()->SetBinLabel(3, "isrDefLo isr:muRfac=2.0");
  PSWeights->GetXaxis()->SetBinLabel(4, "fsrDefLo fsr:muRfac=2.0");

  PDFWeights = fs->make<TH1D>("PDFWeights", "NNPDF31_nnlo_hessian_pdfas (3060000)",103,0,103);
}


ttbbLepJetsAnalyzer::~ttbbLepJetsAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete b_PUWeight;
  delete b_PDFWeight;
  delete b_ScaleWeight;
  delete b_PSWeight;

  delete b_GenConeCatID;
  delete b_GenCone_gJet_pt;
  delete b_GenCone_gJet_eta;
  delete b_GenCone_gJet_phi;
  delete b_GenCone_gJet_e;
  delete b_GenCone_gJetFlavW;

  delete b_Lepton_SF;

  delete b_Jet_pt;
  delete b_Jet_eta;
  delete b_Jet_phi;
  delete b_Jet_e;
  delete b_Jet_Index;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_JES_Up;
  delete b_Jet_JES_Down;
  delete b_Jet_JER_Up;
  delete b_Jet_JER_Nom;
  delete b_Jet_JER_Down;

  delete b_Jet_SF_deepCSV_30;
  delete b_Jet_SF_deepJet_30;

}

//
// member functions
//

// ------------ method called for each event  ------------
void ttbbLepJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_PUWeight     ->clear();
  b_ScaleWeight  ->clear();
  b_PDFWeight    ->clear();
  b_PSWeight     ->clear();

  b_GenConeCatID    ->clear();
  b_GenCone_gJet_pt ->clear();
  b_GenCone_gJet_eta->clear();
  b_GenCone_gJet_phi->clear();
  b_GenCone_gJet_e  ->clear();
  b_GenCone_gJetFlavW->clear();

  b_Lepton_SF->clear();

  b_Jet_pt    ->clear();
  b_Jet_eta   ->clear();
  b_Jet_phi   ->clear();
  b_Jet_e     ->clear();
  b_Jet_Index->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_JES_Up  ->clear();
  b_Jet_JES_Down->clear();
  b_Jet_JER_Up  ->clear();
  b_Jet_JER_Nom ->clear();
  b_Jet_JER_Down->clear();

  b_Jet_deepCSV      ->clear();
  b_Jet_deepJet      ->clear();
  b_Jet_SF_deepCSV_30->clear();
  b_Jet_SF_deepJet_30->clear();
  b_Jet_deepCvsL     ->clear();
  b_Jet_deepCvsB     ->clear();
  b_Jet_deepJetCvsL  ->clear();
  b_Jet_deepJetCvsB  ->clear();

  b_addbjet1_pt  = -1.0;
  b_addbjet1_eta = -1.0;
  b_addbjet1_phi = -1.0;
  b_addbjet1_e   = -1.0;

  b_addbjet2_pt  = -1.0;
  b_addbjet2_eta = -1.0;
  b_addbjet2_phi = -1.0;
  b_addbjet2_e   = -1.0;

  b_mindRjet1_pt  = -1.0;
  b_mindRjet1_eta = -1.0;
  b_mindRjet1_phi = -1.0;
  b_mindRjet1_e   = -1.0;

  b_mindRjet2_pt  = -1.0;
  b_mindRjet2_eta = -1.0;
  b_mindRjet2_phi = -1.0;
  b_mindRjet2_e   = -1.0;

  b_mindR = 999;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  isMC_ = !iEvent.isRealData();

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();

  b_Lepton_relIso = 999;
  b_Lepton_isIso = false;

  EventInfo->Fill(0.5, 1.0);         // Number of Events

  bool METfiltered = false;
  if(isMC_) {

    edm::Handle<int> nTrueVertHandle;
    iEvent.getByToken(nTrueVertToken_, nTrueVertHandle);
    b_nTruePV = *nTrueVertHandle;

    //---------------------------------------------------------------------------
    // PU Info
    //---------------------------------------------------------------------------

    edm::Handle<float> PUWeight;
    iEvent.getByToken(puWeightToken_, PUWeight);
    b_PUWeight->push_back(*PUWeight); // Central

    EventInfo->Fill(8.5, *PUWeight); // Sum of PUWeights

    edm::Handle<float> PUWeight_Up;
    iEvent.getByToken(puUpWeightToken_, PUWeight_Up);
    b_PUWeight->push_back(*PUWeight_Up); //Syst. Up

    edm::Handle<float> PUWeight_Down;
    iEvent.getByToken(puDownWeightToken_, PUWeight_Down);
    b_PUWeight->push_back(*PUWeight_Down); //Syst. Down 

    //---------------------------------------------------------------------------
    // Weights at Generation Level: aMC@NLO
    //---------------------------------------------------------------------------

    edm::Handle<float> genWeight;

    iEvent.getByToken(genWeightToken_, genWeight);
    b_GenWeight = *genWeight;

    EventInfo->Fill(1.5, b_GenWeight); // Sum of aMC@NLO Weights

    //---------------------------------------------------------------------------
    // MET optional filters
    //---------------------------------------------------------------------------
    edm::Handle<int> recoFiltersHandle;
    iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
    METfiltered = *recoFiltersHandle == 0 ? true : false;
  }

  else{
    b_PUWeight  ->push_back(1.0);
    b_GenWeight = 1.0;
  }

  //---------------------------------------------------------------------------
  // Weights for Syst. Scale and PDF: ttbar
  //---------------------------------------------------------------------------
  if(TTbarMC_ == 1 ) {
    // Scale weights
    edm::Handle<std::vector<float>> scaleUpWeightsHandle, scaleDownWeightsHandle;
    iEvent.getByToken(scaleUpWeightToken_,   scaleUpWeightsHandle);
    iEvent.getByToken(scaleDownWeightToken_, scaleDownWeightsHandle);

    b_ScaleWeight->push_back(weight_valid(scaleUpWeightsHandle  ->at(0))); // muR=Nom  muF=Up
    b_ScaleWeight->push_back(weight_valid(scaleDownWeightsHandle->at(0))); // muR=Nom  muF=Down
    b_ScaleWeight->push_back(weight_valid(scaleUpWeightsHandle  ->at(1))); // muR=Up   muF=Nom
    b_ScaleWeight->push_back(weight_valid(scaleUpWeightsHandle  ->at(2))); // muR=Up   muF=Up
    b_ScaleWeight->push_back(weight_valid(scaleDownWeightsHandle->at(1))); // muR=Down muF=Nom
    b_ScaleWeight->push_back(weight_valid(scaleDownWeightsHandle->at(2))); // muR=Down muF=Down
    
    for(unsigned int iscale = 0; iscale< b_ScaleWeight->size(); iscale++)
      ScaleWeights->Fill(iscale, b_ScaleWeight->at(iscale)); 

    // PDF weight : NEED TO BE UPDATE
    edm::Handle<std::vector<float>> PDFWeightsHandle;
    iEvent.getByToken(pdfWeightToken_,   PDFWeightsHandle);

    std::vector<float> tmp_PDFWeight;
    for( auto& w : *PDFWeightsHandle ) tmp_PDFWeight.push_back(w);
    // 9~112: NNPDF31_nnlo_hessian_pdfas
    for( unsigned int i = 9; i < 112; ++i) b_PDFWeight->push_back(tmp_PDFWeight.at(i));
    for( unsigned int i = 0; i < b_PDFWeight->size(); ++i) PDFWeights->Fill(i, b_PDFWeight->at(i));

    // PS weight
    edm::Handle<std::vector<float>> PSWeightsHandle;
    iEvent.getByToken(psWeightToken_, PSWeightsHandle);

    b_PSWeight->push_back(weight_valid(PSWeightsHandle->at(4))); // isrDefHigh
    b_PSWeight->push_back(weight_valid(PSWeightsHandle->at(5))); // fsrDefHigh
    b_PSWeight->push_back(weight_valid(PSWeightsHandle->at(6))); // isrDefLow
    b_PSWeight->push_back(weight_valid(PSWeightsHandle->at(7))); // fsrDefLow

    for(unsigned int iscale = 0; iscale< b_PSWeight->size(); iscale++)
      PSWeights->Fill(iscale, b_PSWeight->at(iscale));
  }
  else{
    b_ScaleWeight->push_back(1.0);
    b_PDFWeight  ->push_back(1.0);
    b_PSWeight   ->push_back(1.0);
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

    b_GenCone_NgJetsW = genttbarConeCat->begin()->NWJets();
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()-> Wquarkflav1());
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()-> Wquarkflav2());

    for (int ijGT=0; ijGT<6; ijGT++){
      b_GenCone_gJet_pt  ->push_back(gJetGenCone[ijGT].Pt());
      b_GenCone_gJet_eta ->push_back(gJetGenCone[ijGT].Eta());
      b_GenCone_gJet_phi ->push_back(gJetGenCone[ijGT].Phi());
      b_GenCone_gJet_e   ->push_back(gJetGenCone[ijGT].E());
    }

    // adding additional b jet four-momentum
    b_addbjet1_pt = genttbarConeCat->begin()->addbJets1().Pt();
    b_addbjet1_eta = genttbarConeCat->begin()->addbJets1().Eta();
    b_addbjet1_phi = genttbarConeCat->begin()->addbJets1().Phi();
    b_addbjet1_e = genttbarConeCat->begin()->addbJets1().E();

    b_addbjet2_pt = genttbarConeCat->begin()->addbJets2().Pt();
    b_addbjet2_eta = genttbarConeCat->begin()->addbJets2().Eta();
    b_addbjet2_phi = genttbarConeCat->begin()->addbJets2().Phi();
    b_addbjet2_e = genttbarConeCat->begin()->addbJets2().E();

    // DR 
    b_DRAddJets = genttbarConeCat->begin()->dRaddJets();

    // mindR b jets
    std::vector<math::XYZTLorentzVector> genbjets = genttbarConeCat->begin()->bJets();
    if(genbjets.size() > 0){
      double mindR=999;
      int index1=-1, index2=-1;
      for(unsigned int i=0; i < genbjets.size(); ++i){
        if(genbjets[i].Pt() < 20 || std::abs(genbjets[i].Eta()) > 2.5) continue;
	for(unsigned int j=i+1; j < genbjets.size(); ++j){
	  if(genbjets[j].Pt() < 20 || std::abs(genbjets[j].Eta()) > 2.5) continue;
	  double tmp = reco::deltaR(genbjets[i], genbjets[j]);
	  if(tmp < mindR){
	    mindR = tmp;
	    index1 = i;
	    index2 = j;
	  }
	}
      }

      if(index1 > -1 && index2 > -1){
        b_mindRjet1_pt  = genbjets[index1].Pt();
        b_mindRjet1_eta = genbjets[index1].Eta();
        b_mindRjet1_phi = genbjets[index1].Phi();
        b_mindRjet1_e   = genbjets[index1].E();

        b_mindRjet2_pt  = genbjets[index2].Pt();
        b_mindRjet2_eta = genbjets[index2].Eta();
        b_mindRjet2_phi = genbjets[index2].Phi();
        b_mindRjet2_e   = genbjets[index2].E();

	b_mindR = mindR;
      }
    }

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
    // Requires ttjj events to be categorized
    if(genttbarConeCat->begin()->NaddJets20() > 1) Isttjj = true;

    if      (Isttjj && genttbarConeCat->begin()->NaddbJets20() > 1) Isttbb = true;
    else if (Isttjj && genttbarConeCat->begin()->NaddbJets20() > 0) Isttb  = true;
    else if (Isttjj && genttbarConeCat->begin()->NaddcJets20() > 1) Isttcc = true;
    else if (Isttjj && genttbarConeCat->begin()->NaddJets20()  > 1) IsttLF = true;
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
	  b_GenLepton_pt  = genttbarConeCat->begin()->lepton1().pt();
	  b_GenLepton_eta = genttbarConeCat->begin()->lepton1().eta();
	  b_GenLepton_phi = genttbarConeCat->begin()->lepton1().phi();
	  b_GenLepton_e   = genttbarConeCat->begin()->lepton1().e();
	  b_GenNu_pt      = genttbarConeCat->begin()->nu1().pt();
	  b_GenNu_eta     = genttbarConeCat->begin()->nu1().eta();
	  b_GenNu_phi     = genttbarConeCat->begin()->nu1().phi();
	  b_GenNu_e       = genttbarConeCat->begin()->nu1().e();
	}
	else{
	  b_GenLepton_pt  = genttbarConeCat->begin()->lepton2().pt();
	  b_GenLepton_eta = genttbarConeCat->begin()->lepton2().eta();
	  b_GenLepton_phi = genttbarConeCat->begin()->lepton2().phi();
	  b_GenLepton_e   = genttbarConeCat->begin()->lepton2().e();
	  b_GenNu_pt      = genttbarConeCat->begin()->nu2().pt();
	  b_GenNu_eta     = genttbarConeCat->begin()->nu2().eta();
	  b_GenNu_phi     = genttbarConeCat->begin()->nu2().phi();
	  b_GenNu_e       = genttbarConeCat->begin()->nu2().e();
	}

        gentree->Fill();

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

  if ( !doLooseLepton_ ) {
    for (unsigned int i = 0; i < electrons->size() ; i++) {
      const cat::Electron & electron = electrons->at(i);

      if     ( IsSelectElectron( electron ) ) selectedElectrons.push_back( electron );
      else if( IsVetoElectron( electron ) )   vetoElectrons.push_back( electron ); // does not Include selected electrons
    }
  }
  else if ( !electrons->empty() ) {
    if ( IsSelectElectron(electrons->at(0)) ) selectedElectrons.push_back( electrons->at(0) );
    for (unsigned int i = 1; i < electrons->size() ; i++) {
      const cat::Electron & electron = electrons->at(i);
      if ( IsVetoElectron(electron) ) vetoElectrons.push_back(electron);
    }
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

  if ( !doLooseLepton_ ) {
    for (unsigned int i = 0; i < muons->size() ; i++) {
      const cat::Muon & muon = muons->at(i);

      if( IsSelectMuon( muon) ) selectedMuons.push_back( muon);
      else if( IsVetoMuon( muon) ) vetoMuons.push_back( muon); // does not Include selected muons
    }
  }
  else if ( !muons->empty() ) {
    if ( IsSelectMuon(muons->at(0)) ) selectedMuons.push_back( muons->at(0) );
    for (unsigned int i = 1; i < muons->size() ; i++) {
      const cat::Muon & muon = muons->at(i);
      if ( IsVetoMuon(muon) ) vetoMuons.push_back(muon);
    }
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  TLorentzVector lepton;
  std::vector<float> *lepton_SF;
  lepton_SF = new std::vector<float>;

  int ch_tag  = 999;

  if(selectedMuons.size()     == 1 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 0 &&
     vetoElectrons.size()     == 0){
    lepton.SetPtEtaPhiE(selectedMuons[0].pt(), selectedMuons[0].eta(), selectedMuons[0].phi(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets
    b_Lepton_relIso = selectedMuons[0].relIso(0.4);
    b_Lepton_isIso = (b_Lepton_relIso < 0.15);

    if(isMC_) {
      // Lepton SF (ID/ISO/Trigger)
      lepton_SF->push_back( SF_muonId_ ( selectedMuons[0].pt(), selectedMuons[0].eta()       )); // [0]-> IdSF
      lepton_SF->push_back( SF_muonId_ ( selectedMuons[0].pt(), selectedMuons[0].eta(),  1.0 )); // [1]-> IdSF+Error
      lepton_SF->push_back( SF_muonId_ ( selectedMuons[0].pt(), selectedMuons[0].eta(), -1.0 )); // [2]-> IdSF-Error
      lepton_SF->push_back( SF_muonIso_( selectedMuons[0].pt(), selectedMuons[0].eta()       )); // [3]-> IsoSF
      lepton_SF->push_back( SF_muonIso_( selectedMuons[0].pt(), selectedMuons[0].eta(),  1.0 )); // [4]-> IsoSF+Error
      lepton_SF->push_back( SF_muonIso_( selectedMuons[0].pt(), selectedMuons[0].eta(), -1.0 )); // [5]-> IsoSF-Error
      lepton_SF->push_back( SF_muonTrg_( selectedMuons[0].pt(), selectedMuons[0].eta()       )); // [6]-> TrgSF
      lepton_SF->push_back( SF_muonTrg_( selectedMuons[0].pt(), selectedMuons[0].eta(),  1.0 )); // [7]-> TrgSF+Error
      lepton_SF->push_back( SF_muonTrg_( selectedMuons[0].pt(), selectedMuons[0].eta(), -1.0 )); // [8]-> TrgSF-Error
    }
  }

  if(selectedMuons.size()     == 0 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 1 &&
     vetoElectrons.size()     == 0){
    ch_tag = 1; //electron + jets
    lepton.SetPtEtaPhiE(selectedElectrons[0].pt(), selectedElectrons[0].eta(), selectedElectrons[0].phi(), selectedElectrons[0].energy());
    b_Lepton_relIso = selectedElectrons.at(0).relIso(0.3);
    if( std::abs(selectedElectrons[0].scEta()) <= 1.479 && b_Lepton_relIso < 0.0287+0.506/static_cast<double>(lepton.Pt())) b_Lepton_isIso = true;
    else if( std::abs(selectedElectrons[0].scEta()) > 1.479 && b_Lepton_relIso < 0.0445+0.963/static_cast<double>(lepton.Pt())) b_Lepton_isIso = true;

    if(isMC_) {
      // Lepton SF (ID/Reco/Zvtx/Trigger)
      lepton_SF->push_back( SF_elecId_  ( selectedElectrons[0].pt(), selectedElectrons[0].scEta()       )); // [0]- > IdSF
      lepton_SF->push_back( SF_elecId_  ( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 )); // [1] -> IdSF+Error
      lepton_SF->push_back( SF_elecId_  ( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 )); // [2] -> IdSF-Error
      lepton_SF->push_back( SF_elecReco_( selectedElectrons[0].pt(), selectedElectrons[0].scEta()       )); // [3] -> RecoSF
      lepton_SF->push_back( SF_elecReco_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 )); // [4] -> RecoSF+Error
      lepton_SF->push_back( SF_elecReco_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 )); // [5] -> RecoSF-Error
      lepton_SF->push_back( SF_elecZvtx_( selectedElectrons[0].pt(), selectedElectrons[0].scEta()       )); // [6] -> ZvtxSF
      lepton_SF->push_back( SF_elecZvtx_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 )); // [7] -> ZvtxSF+Error
      lepton_SF->push_back( SF_elecZvtx_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 )); // [8] -> ZvtxSF-Error
      lepton_SF->push_back( SF_elecTrg_ ( selectedElectrons[0].pt(), selectedElectrons[0].scEta()       )); // [9] -> TrgSF
      lepton_SF->push_back( SF_elecTrg_ ( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 )); // [10]-> TrgSF+Error
      lepton_SF->push_back( SF_elecTrg_ ( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 )); // [11]-> TrgSF-Error
    }
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // HLTrigger
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  bool EvTrigger     = false;
  bool IsTriggerMu   = false;
  bool IsTriggerEl   = false;
  bool IsTriggerElHT = false;

  edm::Handle<int> trigMuFiltersHandle;
  iEvent.getByToken(trigMuFiltersToken_, trigMuFiltersHandle);
  IsTriggerMu = *trigMuFiltersHandle == 0 ? false : true;

  edm::Handle<int> trigElFiltersHandle;
  iEvent.getByToken(trigElFiltersToken_, trigElFiltersHandle);
  IsTriggerEl = *trigElFiltersHandle == 0 ? false : true;

  edm::Handle<int> trigElHTFiltersHandle;
  iEvent.getByToken(trigElHTFiltersToken_, trigElHTFiltersHandle);
  IsTriggerElHT = *trigElHTFiltersHandle == 0 ? false : true;

  if( (ch_tag == 0 && IsTriggerMu) ||
      (ch_tag == 1 && (IsTriggerEl || IsTriggerElHT)) )
    EvTrigger = true;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Fill Tree with events that have ONLY one lepton
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Check Gen Level for ttbar sample
  if (TTbarMC_ >0){
    // Signal ttbar event
    if(TTbarMC_ == 1 && nGenLep != 1) ch_tag = 999;
    // Background ttbar event
    if(TTbarMC_ == 2 && nGenLep == 1) ch_tag = 999;
    // ttbar Category
    if(!IsCat) ch_tag = 999;
  } // if(TTbarMC_ >0)

  if (ch_tag<2 && EvTrigger && !METfiltered){ // Single lepton event

    b_Channel  = ch_tag;

    b_Lepton_pt  = lepton.Pt();
    b_Lepton_eta = lepton.Eta();
    b_Lepton_phi = lepton.Phi();
    b_Lepton_e   = lepton.E();

    b_Lepton_SF  = lepton_SF;

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

    int N_GoodJets = 0;

    // Initialize SF_btag
    // Jet_SF_CSV[Scenario][SystVariations];
    float Jet_SF_deepCSV[1][19];
    for (unsigned int ipTj=0; ipTj<1; ipTj++){
       for (unsigned int iu=0; iu<19; iu++) Jet_SF_deepCSV[ipTj][iu] = 1.0;
    }

    // Run again over all Jets (CSV order)
    for (unsigned int i = 0; i < JetIndex.size() ; i++) {
      const cat::Jet & jet = jets->at(JetIndex[i]);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection (pT>20GeV to take into account SYST Variations)
      if( std::abs(jet.eta()) < 2.4 && jet.pt() > 20. && jet.tightLepVetoJetID() ) goodJet = true;
      // Jet Cleaning
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if( dr_LepJet > 0.4 ) cleanJet = true;

      if(goodJet && cleanJet){
        // Basic variables
        b_Jet_pt ->push_back(jet.pt());
        b_Jet_eta->push_back(jet.eta());
        b_Jet_phi->push_back(jet.phi());
        b_Jet_e  ->push_back(jet.energy());
        b_Jet_Index ->push_back(JetIndex[i]);
	
	// Number of Jets (easy cross check)
        if( jet.pt() > 30. ) N_GoodJets ++;

        // Parton Flavour
        b_Jet_partonFlavour->push_back(jet.partonFlavour());
        b_Jet_hadronFlavour->push_back(jet.hadronFlavour());

        // b-tag discriminant
        float jet_btagDis_deepCSV = 
	  jet.bDiscriminator(BTAG_DeepCSVb) + jet.bDiscriminator(BTAG_DeepCSVbb);
	float jet_btagDis_deepJet = 
	  jet.bDiscriminator(BTAG_DeepJetb)  + 
	  jet.bDiscriminator(BTAG_DeepJetbb) + jet.bDiscriminator(BTAG_DeepJetlepb);
        b_Jet_deepCSV->push_back(jet_btagDis_deepCSV);
        b_Jet_deepJet->push_back(jet_btagDis_deepJet);

        // c-tag discriminant
        float jet_btagDis_deepCvsL = jet.bDiscriminator(BTAG_DeepCSVc)/
	    (jet.bDiscriminator(BTAG_DeepCSVc) + jet.bDiscriminator(BTAG_DeepCSVudsg));
        float jet_btagDis_deepCvsB = jet.bDiscriminator(BTAG_DeepCSVc)/
	    (jet.bDiscriminator(BTAG_DeepCSVc) +
	     jet.bDiscriminator(BTAG_DeepCSVb) + jet.bDiscriminator(BTAG_DeepCSVbb));
	float jet_btagDis_deepJetCvsL = jet.bDiscriminator(BTAG_DeepJetc)/
	    (jet.bDiscriminator(BTAG_DeepJetc)   +
	     jet.bDiscriminator(BTAG_DeepJetuds) + jet.bDiscriminator(BTAG_DeepJetg));
	float jet_btagDis_deepJetCvsB = jet.bDiscriminator(BTAG_DeepJetc)/
	    (jet.bDiscriminator(BTAG_DeepJetc)  + jet.bDiscriminator(BTAG_DeepJetb) +
	     jet.bDiscriminator(BTAG_DeepJetbb) + jet.bDiscriminator(BTAG_DeepJetlepb));
        b_Jet_deepCvsL->push_back(jet_btagDis_deepCvsL);
	b_Jet_deepCvsB->push_back(jet_btagDis_deepCvsB);
	b_Jet_deepJetCvsL->push_back(jet_btagDis_deepJetCvsL);
	b_Jet_deepJetCvsB->push_back(jet_btagDis_deepJetCvsB);

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
          // Saving the central SF and the 18 syst. unc. for:
          if(jet.pt() > 30.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_deepCSV[0][iu] *= SF_deepCSV_.getSF(jet, iu);

        } // if(isMC_)	
      }// if(GoodJets)
    }// for(AllJets)
    
    // Number of Jets (easy cross check)
    b_Jet_Number = N_GoodJets;

    for (unsigned int iu=0; iu<19; iu++) b_Jet_SF_deepCSV_30->push_back(1.0);
    (*b_Jet_SF_deepCSV_30)[0] = Jet_SF_deepCSV[0][0]; //Central
    // To save only the error
    for (unsigned int iu=1; iu<19; iu++){
      (*b_Jet_SF_deepCSV_30)[iu] = std::abs(Jet_SF_deepCSV[0][iu] - Jet_SF_deepCSV[0][0]) ; // Syst. Unc.
    }

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Fill Tree with event at 1 lepton cut level
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    
    tree->Fill();
        
  }// if(ch_tag)
  
}

//------------- Good Muon Selection -----------------------
bool ttbbLepJetsAnalyzer::IsSelectMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Tight selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoTight));

  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 30);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  if ( !doLooseLepton_ ) GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.15 );

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
  GoodElectron &= (i_electron_candidate.pt() > 30);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // Electron cut based selection
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  if ( !doLooseLepton_ ) GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-tight") > 0.0;
  else                   GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-tight-noiso") > 0.0;

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

  // Isolation is already included in the cut-based cuts, it is not needed 
  // if ( !doLooseLepton_ ) {
  //   if ( std::abs(i_electron_candidate.scEta()) <= 1.479)   GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.0695 );
  //   else GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.0821 );
  // }

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
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-veto") > 0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // if ( std::abs(i_electron_candidate.scEta()) <= 1.479) GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.175 );
  // else GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.159 );
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}

float ttbbLepJetsAnalyzer::weight_valid(float weight){
  if( std::isnan(weight) or std::isinf(weight) ) return 1.0;
  if( weight < 0.01 or weight > 100 ) return 1.0;
  return weight;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbbLepJetsAnalyzer);
