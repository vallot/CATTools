// -*- C++ -*-
/**\class fcncLepJetsAnalyzer fcncLepJetsAnalyzer.cc
*/
//
// Package:    fcncLepJetsAnalyzer
// Class:      fcncLepJetsAnalyzer
//
// Original Author: Javier Brochero
// Author:  Jiwon Park
// Created:  2018

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

using namespace cat;

class fcncLepJetsAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit fcncLepJetsAnalyzer(const edm::ParameterSet&);
  ~fcncLepJetsAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  //----------------------------------------------------------------
  bool IsSelectMuon    (const cat::Muon     & i_muon_candidate);
  bool IsVetoMuon      (const cat::Muon     & i_muon_candidate);
  bool IsSelectElectron(const cat::Electron & i_electron_candidate);
  bool IsVetoElectron  (const cat::Electron & i_electron_candidate);
  float weight_valid   (float weight);
  float jer_valid   (float weight);

  bool isMC_;
  bool doLooseLepton_;

  int TTbarMC_; // 0->No ttbar, 1->ttbar Signal
  int TTbarCatMC_;

  // Event Weights
  edm::EDGetTokenT<float>                   genWeightToken_;
  edm::EDGetTokenT<std::vector<float>>      scaleUpWeightToken_, scaleDownWeightToken_;
  edm::EDGetTokenT<std::vector<float>>      psWeightToken_, pdfWeightToken_;
  edm::EDGetTokenT<float>                   puWeightToken_, puUpWeightToken_, puDownWeightToken_;
  edm::EDGetTokenT<double>                  prefweightToken_, prefweightupToken_, prefweightdownToken_;
  // Object Collections
  edm::EDGetTokenT<int>                     genttbarHiggsCatToken_;
  edm::EDGetTokenT<cat::GenTopCollection>   genttbarCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<int>                     pvToken_, nTrueVertToken_, recoFiltersToken_;
  edm::EDGetTokenT<int>                     trigMuFiltersToken_, trigElFiltersToken_, trigElHTFiltersToken_;
// ----------member data ---------------------------

  TTree *tree;
  //TTree *gentree;

  // Event info
  int b_Event, b_Run, b_Lumi_Number;
  float b_GenWeight;
  std::vector<float> *b_ScaleWeight, *b_PDFWeight, *b_PSWeight;

  // PU/Vertices
  std::vector<float> *b_PUWeight;
  int b_nGoodPV, b_nTruePV;
  std::vector<double> *b_PrefireWeight;

  // Channel and Categorization
  int b_Channel, b_GenCone_NgJetsW, b_GenHiggsCatID;
  std::vector<float> *b_GenCone_gJet_pt, *b_GenCone_gJet_eta, *b_GenCone_gJet_phi, *b_GenCone_gJet_e;
  std::vector<int>   *b_GenConeCatID, *b_GenCone_gJetFlavW;

  // GEN Lepton and neutrino
  float b_GenLepton_pt, b_GenLepton_eta, b_GenLepton_phi, b_GenLepton_e;
  float b_GenLepton2_pt, b_GenLepton2_eta, b_GenLepton2_phi, b_GenLepton2_e;
  float b_GenNu_pt, b_GenNu_eta, b_GenNu_phi, b_GenNu_e;
  float b_GenTop1_pt, b_GenTop2_pt, b_TopPtWeight;

  // Leptons and MET
  float b_Lepton_pt, b_Lepton_eta, b_Lepton_phi, b_Lepton_e, b_Lepton_relIso;
  //float b_Lepton_scEta = -99.0;
  float b_MET, b_MET_phi;
  bool  b_Lepton_isIso;
  std::vector<float> *b_Lepton_SF;
  std::vector<float> *b_Electron_Scale;

  //additional jets 
  float b_addbjet1_pt, b_addbjet1_eta, b_addbjet1_phi, b_addbjet1_e; 
  float b_addbjet2_pt, b_addbjet2_eta, b_addbjet2_phi, b_addbjet2_e; 
  float b_Hbjet1_pt, b_Hbjet1_eta, b_Hbjet1_phi, b_Hbjet1_e;
  float b_Hbjet2_pt, b_Hbjet2_eta, b_Hbjet2_phi, b_Hbjet2_e;
  float b_Hbquarkjet1_pt, b_Hbquarkjet1_eta, b_Hbquarkjet1_phi, b_Hbquarkjet1_e;
  float b_Hbquarkjet2_pt, b_Hbquarkjet2_eta, b_Hbquarkjet2_phi, b_Hbquarkjet2_e;
  float b_dRHbb;

  // Jets
  int b_Jet_NJet, b_Jet_NBJetM;
  std::vector<float> *b_Jet_pt, *b_Jet_eta, *b_Jet_phi, *b_Jet_e;
  std::vector<int>   *b_Jet_Index;
  float               b_DRAddJets;

  // Jet Flavour
  std::vector<int> *b_Jet_partonFlavour;
  std::vector<int> *b_Jet_hadronFlavour;

  // JES and JER
  std::vector<float> *b_Jet_JES_Up, *b_Jet_JES_Down;
  std::vector<float> *b_Jet_JER_Up, *b_Jet_JER_Nom, *b_Jet_JER_Down;

  // b-Jet discriminant
  std::vector<float> *b_Jet_deepCSV, *b_Jet_deepJet;
  std::vector<float> *b_Jet_SF_deepCSV_30, *b_Jet_SF_deepJet_30;

  // c-Jet discriminant
  std::vector<float> *b_Jet_deepCvsL, *b_Jet_deepCvsB;
  std::vector<float> *b_Jet_deepJetCvsL, *b_Jet_deepJetCvsB;

  // Histograms: Number of Events and Weights
  TH1D *EventInfo, *ScaleWeights, *PDFWeights, *PSWeights, *TopPtWeight;

  // Scale factor evaluators
  BTagWeightEvaluator SF_deepCSV_;
  ScaleFactorEvaluator SF_muonId_, SF_muonIso_, SF_muonTrg_, SF_elecId_, SF_elecReco_, SF_elecZvtx_, SF_elecTrg_;

};

fcncLepJetsAnalyzer::fcncLepJetsAnalyzer(const edm::ParameterSet& iConfig):
  doLooseLepton_(iConfig.getUntrackedParameter<bool>("doLooseLepton", false)),
  TTbarMC_    (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0)),
  TTbarCatMC_ (iConfig.getUntrackedParameter<int>("TTbarCatLabel", 0))
{
  const auto elecIdSFSet = iConfig.getParameter<edm::ParameterSet>("elecIdSF");
  SF_elecId_.set(elecIdSFSet.getParameter<std::vector<double>>("pt_bins"),
                 elecIdSFSet.getParameter<std::vector<double>>("eta_bins"),
                 elecIdSFSet.getParameter<std::vector<double>>("values"),
                 elecIdSFSet.getParameter<std::vector<double>>("errors"));
  const auto elecRecoSFSet = iConfig.getParameter<edm::ParameterSet>("elecRecoSF");
  SF_elecReco_.set(elecRecoSFSet.getParameter<std::vector<double>>("pt_bins"),
                   elecRecoSFSet.getParameter<std::vector<double>>("eta_bins"),
                   elecRecoSFSet.getParameter<std::vector<double>>("values"),
                   elecRecoSFSet.getParameter<std::vector<double>>("errors"));
  const auto elecZvtxSFSet = iConfig.getParameter<edm::ParameterSet>("elecZvtxSF");
  SF_elecZvtx_.set(elecZvtxSFSet.getParameter<std::vector<double>>("pt_bins"),
                   elecZvtxSFSet.getParameter<std::vector<double>>("eta_bins"),
                   elecZvtxSFSet.getParameter<std::vector<double>>("values"),
                   elecZvtxSFSet.getParameter<std::vector<double>>("errors"));
  const auto elecTrgSFSet = iConfig.getParameter<edm::ParameterSet>("elecTrgSF");
  SF_elecTrg_.set(elecTrgSFSet.getParameter<std::vector<double>>("pt_bins"),
                  elecTrgSFSet.getParameter<std::vector<double>>("eta_bins"),
                  elecTrgSFSet.getParameter<std::vector<double>>("values"),
                  elecTrgSFSet.getParameter<std::vector<double>>("errors"));
  const auto muonIdSFSet = iConfig.getParameter<edm::ParameterSet>("muonIdSF");
  SF_muonId_.set(muonIdSFSet.getParameter<std::vector<double>>("pt_bins"),
                 muonIdSFSet.getParameter<std::vector<double>>("abseta_bins"),
                 muonIdSFSet.getParameter<std::vector<double>>("values"),
                 muonIdSFSet.getParameter<std::vector<double>>("errors"));
  const auto muonIsoSFSet = iConfig.getParameter<edm::ParameterSet>("muonIsoSF");
  SF_muonIso_.set(muonIsoSFSet.getParameter<std::vector<double>>("pt_bins"),
                  muonIsoSFSet.getParameter<std::vector<double>>("abseta_bins"),
                  muonIsoSFSet.getParameter<std::vector<double>>("values"),
                  muonIsoSFSet.getParameter<std::vector<double>>("errors"));
  const auto muonTrgSFSet = iConfig.getParameter<edm::ParameterSet>("muonTrgSF");
  SF_muonTrg_.set(muonTrgSFSet.getParameter<std::vector<double>>("pt_bins"),
                  muonTrgSFSet.getParameter<std::vector<double>>("abseta_bins"),
                  muonTrgSFSet.getParameter<std::vector<double>>("values"),
                  muonTrgSFSet.getParameter<std::vector<double>>("errors"));

  SF_deepCSV_.initCSVWeight(false, "deepcsv");

  // Weights
  auto genWeightLabel = iConfig.getParameter<edm::InputTag>("genWeightLabel");
  genWeightToken_       = consumes<float>             (edm::InputTag(genWeightLabel.label()));
  // PDF, Scale, PS
  pdfWeightToken_       = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "pdf"));
  scaleUpWeightToken_   = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaleup"));
  scaleDownWeightToken_ = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaledown"));
  psWeightToken_        = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "ps"));
  // PileUp
  auto puWeightLabel    = iConfig.getParameter<edm::InputTag>("puWeightLabel");
  puWeightToken_        = consumes<float>(edm::InputTag(puWeightLabel.label()));
  puUpWeightToken_      = consumes<float>(edm::InputTag(puWeightLabel.label(),"up"));
  puDownWeightToken_    = consumes<float>(edm::InputTag(puWeightLabel.label(),"dn"));
  // Prefire (16 and 17 only)
  prefweightToken_      = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"));
  prefweightupToken_    = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
  prefweightdownToken_  = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
  // GEN and ttbar Categorization
  genttbarCatToken_      = consumes<cat::GenTopCollection>(iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));
  genttbarHiggsCatToken_ = consumes<int>                  (iConfig.getParameter<edm::InputTag>("genHiggsCatLabel"));
  // Object Collections
  muonToken_        = consumes<cat::MuonCollection>    (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_    = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_         = consumes<cat::JetCollection>     (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_         = consumes<cat::METCollection>     (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_          = consumes<int>                    (iConfig.getParameter<edm::InputTag>("pvLabel"));
  // Trigger
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  trigMuFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMuFilters"));
  trigElFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigElFilters"));
  trigElHTFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigElHTFilters"));

  // PU = 0 prescription
  nTrueVertToken_    = consumes<int>(iConfig.getParameter<edm::InputTag>("nTrueVertLabel"));
 
  b_PUWeight       = new std::vector<float>;
  b_PrefireWeight  = new std::vector<double>;
  b_PDFWeight      = new std::vector<float>;
  b_ScaleWeight    = new std::vector<float>;
  b_PSWeight       = new std::vector<float>;
  b_Lepton_SF      = new std::vector<float>;
  b_Electron_Scale = new std::vector<float>;

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
  b_Jet_deepJetCvsL      = new std::vector<float>;
  b_Jet_deepJetCvsB      = new std::vector<float>;
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
  tree->Branch("prefireweight", "std::vector<double>",&b_PrefireWeight);
  tree->Branch("pdfweight",     "std::vector<float>", &b_PDFWeight );
  tree->Branch("scaleweight",   "std::vector<float>", &b_ScaleWeight );
  tree->Branch("psweight",      "std::vector<float>", &b_PSWeight );
  tree->Branch("topptweight",   &b_TopPtWeight,       "topptweight/F");

  tree->Branch("MET",           &b_MET,        "MET/F");
  tree->Branch("MET_phi",       &b_MET_phi,    "MET_phi/F");
  tree->Branch("lepton_pt",     &b_Lepton_pt,  "lepton_pt/F");
  tree->Branch("lepton_eta",    &b_Lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi",    &b_Lepton_phi, "lepton_phi/F");
  tree->Branch("lepton_e" ,     &b_Lepton_e,   "lepton_e/F" );
  tree->Branch("lepton_SF",     "std::vector<float>", &b_Lepton_SF );
  tree->Branch("lepton_scale",  "std::vector<float>", &b_Electron_Scale );
  tree->Branch("lepton_relIso", &b_Lepton_relIso, "lepton_relIso/F");
  tree->Branch("lepton_isIso",  &b_Lepton_isIso,  "lepton_isIso/O");
  //tree->Branch("lepton_scEta",  &b_Lepton_scEta,  "lepton_scEta/F");

  tree->Branch("jet_pt",            "std::vector<float>", &b_Jet_pt);
  tree->Branch("jet_eta",           "std::vector<float>", &b_Jet_eta);
  tree->Branch("jet_phi",           "std::vector<float>", &b_Jet_phi);
  tree->Branch("jet_e" ,            "std::vector<float>", &b_Jet_e );
  tree->Branch("jet_index" ,        "std::vector<int>",   &b_Jet_Index );
  tree->Branch("jet_deepCSV",       "std::vector<float>", &b_Jet_deepCSV );
  tree->Branch("jet_deepJet",       "std::vector<float>", &b_Jet_deepJet );
  tree->Branch("jet_SF_deepCSV_30", "std::vector<float>", &b_Jet_SF_deepCSV_30 );
  tree->Branch("jet_SF_deepJet_30", "std::vector<float>", &b_Jet_SF_deepJet_30 );
  tree->Branch("jet_deepCvsL",      "std::vector<float>", &b_Jet_deepCvsL );
  tree->Branch("jet_deepCvsB",      "std::vector<float>", &b_Jet_deepCvsB );
  tree->Branch("jet_deepJetCvsL",   "std::vector<float>", &b_Jet_deepJetCvsL );
  tree->Branch("jet_deepJetCvsB",   "std::vector<float>", &b_Jet_deepJetCvsB );

  tree->Branch("jet_njet",          &b_Jet_NJet,          "jet_njet/I" );
  tree->Branch("jet_nbjetm",        &b_Jet_NBJetM,        "jet_nbjetm/I" );
  tree->Branch("jet_partonFlavour", "std::vector<int>",   &b_Jet_partonFlavour );
  tree->Branch("jet_hadronFlavour", "std::vector<int>",   &b_Jet_hadronFlavour );
  tree->Branch("jet_JES_Up",        "std::vector<float>", &b_Jet_JES_Up );
  tree->Branch("jet_JES_Down",      "std::vector<float>", &b_Jet_JES_Down );
  tree->Branch("jet_JER_Up",        "std::vector<float>", &b_Jet_JER_Up );
  tree->Branch("jet_JER_Nom",       "std::vector<float>", &b_Jet_JER_Nom );
  tree->Branch("jet_JER_Down",      "std::vector<float>", &b_Jet_JER_Down );

  //add bjets from Higgs
  tree->Branch("Hbjet1_pt",  &b_Hbjet1_pt,  "Hbjet1_pt/F");
  tree->Branch("Hbjet1_eta", &b_Hbjet1_eta, "Hbjet1_eta/F");
  tree->Branch("Hbjet1_phi", &b_Hbjet1_phi, "Hbjet1_phi/F");
  tree->Branch("Hbjet1_e",   &b_Hbjet1_e,   "Hbjet1_e/F");
  tree->Branch("Hbjet2_pt",  &b_Hbjet2_pt,  "Hbjet2_pt/F");
  tree->Branch("Hbjet2_eta", &b_Hbjet2_eta, "Hbjet2_eta/F");
  tree->Branch("Hbjet2_phi", &b_Hbjet2_phi, "Hbjet2_phi/F");
  tree->Branch("Hbjet2_e",   &b_Hbjet2_e,   "Hbjet2_e/F");
  tree->Branch("dRHbb",      &b_dRHbb,      "dRHbb/F");

  tree->Branch("Hbquarkjet1_pt",  &b_Hbquarkjet1_pt,  "Hbquarkjet1_pt/F");
  tree->Branch("Hbquarkjet1_eta", &b_Hbquarkjet1_eta, "Hbquarkjet1_eta/F");
  tree->Branch("Hbquarkjet1_phi", &b_Hbquarkjet1_phi, "Hbquarkjet1_phi/F");
  tree->Branch("Hbquarkjet1_e",   &b_Hbquarkjet1_e,   "Hbquarkjet1_e/F");
  tree->Branch("Hbquarkjet2_pt",  &b_Hbquarkjet2_pt,  "Hbquarkjet2_pt/F");
  tree->Branch("Hbquarkjet2_eta", &b_Hbquarkjet2_eta, "Hbquarkjet2_eta/F");
  tree->Branch("Hbquarkjet2_phi", &b_Hbquarkjet2_phi, "Hbquarkjet2_phi/F");
  tree->Branch("Hbquarkjet2_e",   &b_Hbquarkjet2_e,   "Hbquarkjet2_e/F");

  // GEN Variables (only ttbarSignal)
  if(TTbarMC_ == 1){
 
    tree->Branch("genconecatid",      "std::vector<int>",   &b_GenConeCatID);
    tree->Branch("gencone_gjet_pt",   "std::vector<float>", &b_GenCone_gJet_pt);
    tree->Branch("gencone_gjet_eta",  "std::vector<float>", &b_GenCone_gJet_eta);
    tree->Branch("gencone_gjet_phi",  "std::vector<float>", &b_GenCone_gJet_phi);
    tree->Branch("gencone_gjet_e",    "std::vector<float>", &b_GenCone_gJet_e);
    tree->Branch("gencone_gJetFlavW", "std::vector<int>",   &b_GenCone_gJetFlavW);
    tree->Branch("gencone_NgjetsW",   &b_GenCone_NgJetsW,   "gencone_NgjetsW/I");
    tree->Branch("draddjets",         &b_DRAddJets,         "draddjets/F");
    tree->Branch("genhiggscatid",     &b_GenHiggsCatID,     "genhiggscatid/I");

    tree->Branch("genlepton_pt",  &b_GenLepton_pt,  "genlepton_pt/F");
    tree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    tree->Branch("genlepton_phi", &b_GenLepton_phi, "genlepton_phi/F");
    tree->Branch("genlepton_e",   &b_GenLepton_e,   "genlepton_e/F");
    tree->Branch("gennu_pt",      &b_GenNu_pt,      "gennu_pt/F");
    tree->Branch("gennu_eta",     &b_GenNu_eta,     "gennu_eta/F");
    tree->Branch("gennu_phi",     &b_GenNu_phi,     "gennu_phi/F");
    tree->Branch("gennu_e",       &b_GenNu_e,       "gennu_e/F");
    tree->Branch("gentop1_pt",    &b_GenTop1_pt,    "gentop1_pt/F");
    tree->Branch("gentop2_pt",    &b_GenTop2_pt,    "gentop2_pt/F");

    tree->Branch("addbjet1_pt",  &b_addbjet1_pt,  "addbjet1_pt/F");
    tree->Branch("addbjet1_eta", &b_addbjet1_eta, "addbjet1_eta/F");
    tree->Branch("addbjet1_phi", &b_addbjet1_phi, "addbjet1_phi/F");
    tree->Branch("addbjet1_e",   &b_addbjet1_e,   "addbjet1_e/F");
    tree->Branch("addbjet2_pt",  &b_addbjet2_pt,  "addbjet2_pt/F");
    tree->Branch("addbjet2_eta", &b_addbjet2_eta, "addbjet2_eta/F");
    tree->Branch("addbjet2_phi", &b_addbjet2_phi, "addbjet2_phi/F");
    tree->Branch("addbjet2_e",   &b_addbjet2_e,   "addbjet2_e/F");

    /*
    //gen tree
    if( TTbarCatMC_ == 1 ){
      gentree = fs->make<TTree>("gentree", "TopGENTree");
      
      gentree->Branch("genweight",     &b_GenWeight,     "genweight/F");
      gentree->Branch("scaleweight",   "std::vector<float>", &b_ScaleWeight );
      gentree->Branch("genconecatid",  "std::vector<int>",   &b_GenConeCatID);
      gentree->Branch("draddjets",     &b_DRAddJets,     "draddjets/F");
      gentree->Branch("genlepton_pt",  &b_GenLepton_pt,  "genlepton_pt/F");
      gentree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
      gentree->Branch("genlepton_phi", &b_GenLepton_phi, "genlepton_phi/F");
      gentree->Branch("genlepton_e",   &b_GenLepton_e,   "genlepton_e/F");
      gentree->Branch("genlepton2_pt", &b_GenLepton2_pt, "genlepton2_pt/F");
      gentree->Branch("genlepton2_eta",&b_GenLepton2_eta,"genlepton2_eta/F");
      gentree->Branch("genlepton2_phi",&b_GenLepton2_phi,"genlepton2_phi/F");
      gentree->Branch("genlepton2_e",  &b_GenLepton2_e,  "genlepton2_e/F");

      gentree->Branch("addbjet1_pt",  &b_addbjet1_pt,  "addbjet1_pt/F");
      gentree->Branch("addbjet1_eta", &b_addbjet1_eta, "addbjet1_eta/F");
      gentree->Branch("addbjet1_phi", &b_addbjet1_phi, "addbjet1_phi/F");
      gentree->Branch("addbjet1_e",   &b_addbjet1_e,   "addbjet1_e/F");
      gentree->Branch("addbjet2_pt",  &b_addbjet2_pt,  "addbjet2_pt/F");
      gentree->Branch("addbjet2_eta", &b_addbjet2_eta, "addbjet2_eta/F");
      gentree->Branch("addbjet2_phi", &b_addbjet2_phi, "addbjet2_phi/F");
      gentree->Branch("addbjet2_e",   &b_addbjet2_e,   "addbjet2_e/F");
    }
    */
  }

  EventInfo = fs->make<TH1D>("EventInfo","Event Information",3,0,3);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
  EventInfo->GetXaxis()->SetBinLabel(3,"Sum of PU Weights");

  ScaleWeights = fs->make<TH1D>("ScaleWeights","Event Weights",6,0,6);
  ScaleWeights->GetXaxis()->SetBinLabel(1,"muR=Nom  muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(2,"muR=Nom  muF=Down");
  ScaleWeights->GetXaxis()->SetBinLabel(3,"muR=Up   muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(4,"muR=Down muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(5,"muR=Up   muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(6,"muR=Down muF=Down");

  PSWeights = fs->make<TH1D>("PSWeights","PS Weights",4,0,4);
  PSWeights->GetXaxis()->SetBinLabel(1,"isrDefHi isr:muRfac=0.5");
  PSWeights->GetXaxis()->SetBinLabel(2,"fsrDefHi fsr:muRfac=0.5");
  PSWeights->GetXaxis()->SetBinLabel(3,"isrDefLo isr:muRfac=2.0");
  PSWeights->GetXaxis()->SetBinLabel(4,"fsrDefLo fsr:muRfac=2.0");

  PDFWeights = fs->make<TH1D>("PDFWeights","PDF4LHC15_nnlo_100_pdfas (91200)",103,0,103);
  TopPtWeight = fs->make<TH1D>("TopPtWeight","Top Pt Reweight",1,0,1);
}


fcncLepJetsAnalyzer::~fcncLepJetsAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete b_PUWeight;
  delete b_PrefireWeight;
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
  delete b_Electron_Scale;
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
  delete b_Jet_deepCSV;
  delete b_Jet_deepJet;
  delete b_Jet_deepCvsL;
  delete b_Jet_deepCvsB;
  delete b_Jet_deepJetCvsL;
  delete b_Jet_deepJetCvsB;

}

// ------------ method called for each event  ------------
void fcncLepJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_PUWeight     ->clear();
  b_PrefireWeight->clear();
  b_ScaleWeight  ->clear();
  b_PSWeight     ->clear();
  b_PDFWeight    ->clear();

  b_GenConeCatID     ->clear();
  b_GenCone_gJet_pt  ->clear();
  b_GenCone_gJet_eta ->clear();
  b_GenCone_gJet_phi ->clear();
  b_GenCone_gJet_e   ->clear();
  b_GenCone_gJetFlavW->clear();

  b_Lepton_SF->clear();
  b_Electron_Scale->clear();

  b_Jet_pt   ->clear();
  b_Jet_eta  ->clear();
  b_Jet_phi  ->clear();
  b_Jet_e    ->clear();
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

  b_addbjet1_pt  = b_addbjet1_e   = b_addbjet2_pt  = b_addbjet2_e   = -1.0 ;
  b_addbjet1_eta = b_addbjet1_phi = b_addbjet2_eta = b_addbjet2_phi = -10.0;
  b_Hbjet1_pt  = b_Hbjet1_e   = b_Hbjet2_pt  = b_Hbjet2_e   = -1.0;
  b_Hbjet1_eta = b_Hbjet1_phi = b_Hbjet2_eta = b_Hbjet2_phi = -10.0;
  b_Hbquarkjet1_pt  = b_Hbquarkjet1_e   = b_Hbquarkjet2_pt  = b_Hbquarkjet2_e   = -1.0;
  b_Hbquarkjet1_eta = b_Hbquarkjet1_phi = b_Hbquarkjet2_eta = b_Hbquarkjet2_phi = -10.0;
  b_dRHbb = -1.0;
  b_TopPtWeight = 1.0;

  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  isMC_ = !iEvent.isRealData();

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();
  b_Lepton_relIso = 999;
  b_Lepton_isIso = false;

  EventInfo->Fill(0.5, 1.0);         // Number of Events

  bool METfiltered = false;
  if( isMC_ ) {

    edm::Handle<int> nTrueVertHandle;
    iEvent.getByToken( nTrueVertToken_, nTrueVertHandle );
    b_nTruePV = *nTrueVertHandle;

    edm::Handle<float> PUWeight;
    iEvent.getByToken(puWeightToken_, PUWeight);
    b_PUWeight->push_back(*PUWeight); // Central
    EventInfo->Fill(2.5, *PUWeight); // Sum of PUWeights

    edm::Handle<float> PUWeight_Up;
    iEvent.getByToken(puUpWeightToken_, PUWeight_Up);
    b_PUWeight->push_back(*PUWeight_Up); //Syst. Up

    edm::Handle<float> PUWeight_Down;
    iEvent.getByToken(puDownWeightToken_, PUWeight_Down);
    b_PUWeight->push_back(*PUWeight_Down); //Syst. Up

    edm::Handle<float> genWeight;
    iEvent.getByToken(genWeightToken_, genWeight);
    b_GenWeight = *genWeight;
    EventInfo->Fill(1.5, b_GenWeight); // Sum of aMC@NLO Weights

    edm::Handle<double> prefweight;
    iEvent.getByToken(prefweightToken_, prefweight );
    b_PrefireWeight->push_back(*prefweight);

    edm::Handle<double> prefweight_up;
    iEvent.getByToken(prefweightupToken_, prefweight_up );
    b_PrefireWeight->push_back(*prefweight_up);

    edm::Handle<double> prefweight_down;
    iEvent.getByToken(prefweightdownToken_, prefweight_down );
    b_PrefireWeight->push_back(*prefweight_down);

    //---------------------------------------------------------------------------
    // MET optional filters 
    //---------------------------------------------------------------------------
    edm::Handle<int> recoFiltersHandle;
    iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
    METfiltered = *recoFiltersHandle == 0 ? true : false;

  }
  else{
    b_PUWeight->push_back(1.0);
    b_ScaleWeight->push_back(1.0);
    b_PDFWeight->push_back(1.0);
    b_PSWeight->push_back(1.0);
    b_GenWeight = 1.0;
    b_nTruePV = 0;
  }

  //---------------------------------------------------------------------------
  // Weights for Syst. Scale and PDF: ttbar
  //---------------------------------------------------------------------------
  if( TTbarMC_ == 1 && TTbarCatMC_ < 4 ) { //0~3 for powheg TT, 4 for ST/TT FCNC

    // muR/muF Scale Weights, exclude nan, empty or crazy valuesa
    edm::Handle<std::vector<float>> scaleUpWeightsHandle, scaleDownWeightsHandle;
    iEvent.getByToken(scaleUpWeightToken_,   scaleUpWeightsHandle);
    iEvent.getByToken(scaleDownWeightToken_, scaleDownWeightsHandle);

    b_ScaleWeight->push_back( weight_valid(scaleUpWeightsHandle  ->at(0)) ); // muR=Nom  muF=Up
    b_ScaleWeight->push_back( weight_valid(scaleDownWeightsHandle->at(0)) ); // muR=Nom  muF=Down
    b_ScaleWeight->push_back( weight_valid(scaleUpWeightsHandle  ->at(1)) ); // muR=Up   muF=Nom
    b_ScaleWeight->push_back( weight_valid(scaleDownWeightsHandle->at(1)) ); // muR=Down muF=Nom
    b_ScaleWeight->push_back( weight_valid(scaleUpWeightsHandle  ->at(2)) ); // muR=Up   muF=Up
    b_ScaleWeight->push_back( weight_valid(scaleDownWeightsHandle->at(2)) ); // muR=Down muF=Down

    for( unsigned int iscale = 0; iscale< b_ScaleWeight->size(); iscale++ )
      ScaleWeights->Fill(iscale, b_ScaleWeight->at(iscale)); 

    // PDF weight
    edm::Handle<std::vector<float>> PDFWeightsHandle;
    iEvent.getByToken(pdfWeightToken_, PDFWeightsHandle);

    std::vector<float> tmp_PDFWeight;
    for( auto& w : *PDFWeightsHandle ) tmp_PDFWeight.push_back(w);
    for( unsigned int i = 569; i < 672; i++)
      b_PDFWeight->push_back(tmp_PDFWeight.at(i));

    for( unsigned int iscale = 0; iscale< b_PDFWeight->size(); iscale++ ) // 578~680: PDF4LHC15_nnlo_100_pdfas, and 0~8 is scale weight
      PDFWeights->Fill(iscale, b_PDFWeight->at(iscale));

    // PS weight
    edm::Handle<std::vector<float>> PSWeightsHandle;
    iEvent.getByToken(psWeightToken_, PSWeightsHandle);

    b_PSWeight->push_back( weight_valid(PSWeightsHandle->at(4)) ); //isrDefHi
    b_PSWeight->push_back( weight_valid(PSWeightsHandle->at(5)) ); //fsrDefHi
    b_PSWeight->push_back( weight_valid(PSWeightsHandle->at(6)) ); //isrDefLo
    b_PSWeight->push_back( weight_valid(PSWeightsHandle->at(7)) ); //fsrDefLo

    for( unsigned int iscale = 0; iscale< b_PSWeight->size(); iscale++ )
      PSWeights->Fill(iscale, b_PSWeight->at(iscale));
  }
  else{
    b_ScaleWeight->push_back(1.0);
    b_PDFWeight->push_back(1.0);
    b_PSWeight->push_back(1.0);
  }

  //---------------------------------------------------------------------------
  // Generated Particles (For Pythia8)
  //---------------------------------------------------------------------------
  bool IsCat = false;
  if( TTbarMC_ > 0 ) {
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

    b_GenConeCatID->push_back(genttbarConeCat->begin()->semiLeptonic(0));
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NJets20());
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NbJets20());
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NcJets20());
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NbJets20NoTop());
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NaddJets20());
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NaddbJets20());
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NaddcJets20());

    //-----------------
    //add b jets from H
    //-----------------
    b_Hbjet1_pt  = genttbarConeCat->begin()->HbJets1().Pt();
    b_Hbjet1_eta = genttbarConeCat->begin()->HbJets1().Eta();
    b_Hbjet1_phi = genttbarConeCat->begin()->HbJets1().Phi();
    b_Hbjet1_e   = genttbarConeCat->begin()->HbJets1().E();
    b_Hbjet2_pt  = genttbarConeCat->begin()->HbJets2().Pt();
    b_Hbjet2_eta = genttbarConeCat->begin()->HbJets2().Eta();
    b_Hbjet2_phi = genttbarConeCat->begin()->HbJets2().Phi();
    b_Hbjet2_e   = genttbarConeCat->begin()->HbJets2().E();
    b_dRHbb      = genttbarConeCat->begin()->dRbJetsFromHiggs();
    b_Hbquarkjet1_pt  = genttbarConeCat->begin()->HbquarkJets1().Pt();
    b_Hbquarkjet1_eta = genttbarConeCat->begin()->HbquarkJets1().Eta();
    b_Hbquarkjet1_phi = genttbarConeCat->begin()->HbquarkJets1().Phi();
    b_Hbquarkjet1_e   = genttbarConeCat->begin()->HbquarkJets1().E();
    b_Hbquarkjet2_pt  = genttbarConeCat->begin()->HbquarkJets2().Pt();
    b_Hbquarkjet2_eta = genttbarConeCat->begin()->HbquarkJets2().Eta();
    b_Hbquarkjet2_phi = genttbarConeCat->begin()->HbquarkJets2().Phi();
    b_Hbquarkjet2_e   = genttbarConeCat->begin()->HbquarkJets2().E();

    // additional Jets Information (p4)
    math::XYZTLorentzVector gJetGenCone[7];
    gJetGenCone[0] = genttbarConeCat->begin()->bJetsFromTop1();
    gJetGenCone[1] = genttbarConeCat->begin()->bJetsFromTop2();
    gJetGenCone[2] = genttbarConeCat->begin()->JetsFromW1();
    gJetGenCone[3] = genttbarConeCat->begin()->JetsFromW2();
    gJetGenCone[4] = genttbarConeCat->begin()->addJets1();
    gJetGenCone[5] = genttbarConeCat->begin()->addJets2();
    gJetGenCone[6] = genttbarConeCat->begin()->FCNCupJet();

    b_GenCone_NgJetsW = genttbarConeCat->begin()->NWJets();
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()->Wquarkflav1());
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()->Wquarkflav2());

    for( int ijGT=0; ijGT<7; ijGT++ ){
      b_GenCone_gJet_pt  ->push_back(gJetGenCone[ijGT].Pt());
      b_GenCone_gJet_eta ->push_back(gJetGenCone[ijGT].Eta());
      b_GenCone_gJet_phi ->push_back(gJetGenCone[ijGT].Phi());
      b_GenCone_gJet_e   ->push_back(gJetGenCone[ijGT].E());
    }

    b_GenTop1_pt = genttbarConeCat->begin()->topquark1().Pt();
    b_GenTop2_pt = genttbarConeCat->begin()->topquark2().Pt();
    float topPtWeight1 = 1.0;
    float topPtWeight2 = 1.0;
    if( b_GenTop1_pt > 0 and b_GenTop1_pt < 800 ){
      topPtWeight1 = TMath::Exp(0.0615-0.0005*b_GenTop1_pt);
    }
    if( b_GenTop2_pt > 0 and b_GenTop2_pt < 800 ){
      topPtWeight1 = TMath::Exp(0.0615-0.0005*b_GenTop2_pt);
    }
    if( TTbarCatMC_ < 4 ){
      b_TopPtWeight = topPtWeight1 * topPtWeight2;
    }

    // adding additional b jet four-momentum
    b_addbjet1_pt  = genttbarConeCat->begin()->addbJets1().Pt();
    b_addbjet1_eta = genttbarConeCat->begin()->addbJets1().Eta();
    b_addbjet1_phi = genttbarConeCat->begin()->addbJets1().Phi();
    b_addbjet1_e   = genttbarConeCat->begin()->addbJets1().E();
    b_addbjet2_pt  = genttbarConeCat->begin()->addbJets2().Pt();
    b_addbjet2_eta = genttbarConeCat->begin()->addbJets2().Eta();
    b_addbjet2_phi = genttbarConeCat->begin()->addbJets2().Phi();
    b_addbjet2_e   = genttbarConeCat->begin()->addbJets2().E();
    b_DRAddJets    = genttbarConeCat->begin()->dRaddJets();

    
    //---------------------------------------------------------------------------
    // Using the GenChannel from GenTop categorization
    //---------------------------------------------------------------------------
    // Category
    int ttbarCAT = b_GenHiggsCatID%100;
    bool Isttbb = false;
    bool Isttcc = false;
    bool IsttLF = false;

    // Categorization based in the Full Ph-Sp
    // Requires ttjj events to be categorized
    if      ( ttbarCAT == 51 or ttbarCAT == 52 or ttbarCAT == 53 or ttbarCAT == 54 or ttbarCAT == 55 ) Isttbb = true;
    else if ( ttbarCAT == 41 or ttbarCAT == 42 or ttbarCAT == 43 or ttbarCAT == 44 or ttbarCAT == 45 ) Isttcc = true;
    else    IsttLF = true;

    if( TTbarCatMC_ == 0 || TTbarCatMC_ == 4 ) IsCat = true;//no cat. for fcnc, ttV
    if( Isttbb && TTbarCatMC_ == 1 ) IsCat = true;
    if( Isttcc && TTbarCatMC_ == 2 ) IsCat = true;
    if( IsttLF && TTbarCatMC_ == 3 ) IsCat = true;

    if( isMC_ && TTbarMC_== 1 && IsCat ){
      if( genttbarConeCat->begin()->lepton1().pt() != 0. ){
        b_GenLepton_pt  = genttbarConeCat->begin()->lepton1().Pt();
        b_GenLepton_eta = genttbarConeCat->begin()->lepton1().Eta();
        b_GenLepton_phi = genttbarConeCat->begin()->lepton1().Phi();
        b_GenLepton_e   = genttbarConeCat->begin()->lepton1().E();
        b_GenNu_pt      = genttbarConeCat->begin()->nu1().Pt();
        b_GenNu_eta     = genttbarConeCat->begin()->nu1().Eta();
        b_GenNu_phi     = genttbarConeCat->begin()->nu1().Phi();
        b_GenNu_e       = genttbarConeCat->begin()->nu1().E();

        if( genttbarConeCat->begin()->lepton2().pt() != 0. ){
          b_GenLepton2_pt = genttbarConeCat->begin()->lepton2().Pt();
          b_GenLepton2_eta= genttbarConeCat->begin()->lepton2().Eta();
          b_GenLepton2_phi= genttbarConeCat->begin()->lepton2().Phi();
          b_GenLepton2_e  = genttbarConeCat->begin()->lepton2().E();
        }
      }
      else{
        b_GenLepton_pt  = genttbarConeCat->begin()->lepton2().Pt();
        b_GenLepton_eta = genttbarConeCat->begin()->lepton2().Eta();
        b_GenLepton_phi = genttbarConeCat->begin()->lepton2().Phi();
        b_GenLepton_e   = genttbarConeCat->begin()->lepton2().E();
        b_GenNu_pt      = genttbarConeCat->begin()->nu2().Pt();
        b_GenNu_eta     = genttbarConeCat->begin()->nu2().Eta();
        b_GenNu_phi     = genttbarConeCat->begin()->nu2().Phi();
        b_GenNu_e       = genttbarConeCat->begin()->nu2().E();
      }

      //if( TTbarCatMC_ == 1 ) gentree->Fill();
    }// if(GENTTbarMCTree_)
  } // if(TTbarMC==0)

  TopPtWeight->Fill(0.5, b_TopPtWeight);

  //---------------------------------------------------------------------------
  // Primary Vertex Info
  //---------------------------------------------------------------------------
  edm::Handle<int> pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );
  b_nGoodPV = *pvHandle;

  //---------------------------------------------------------------------------
  // Missing E_T
  //---------------------------------------------------------------------------
  Handle<cat::METCollection> MET;
  iEvent.getByToken(metToken_, MET);

  b_MET     = MET->begin()->pt();
  b_MET_phi = MET->begin()->phi();

  //---------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------
  cat::ElectronCollection selectedElectrons;
  cat::ElectronCollection vetoElectrons;

  Handle<cat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  if( !doLooseLepton_ ){
    for( unsigned int i = 0; i < electrons->size() ; i++ ){
      const cat::Electron & electron = electrons->at(i);
      if     ( IsSelectElectron( electron ) ) selectedElectrons.push_back( electron );
      else if( IsVetoElectron( electron ) ) vetoElectrons.push_back( electron ); // does not Include selected electrons
    }
  }
  else if( !electrons->empty() ){
    if( IsSelectElectron(electrons->at(0)) ) selectedElectrons.push_back( electrons->at(0) );
    for(unsigned int i = 1; i < electrons->size() ; i++){
      const cat::Electron & electron = electrons->at(i);
      if( IsVetoElectron(electron) ) vetoElectrons.push_back(electron);
    }
  }

  //---------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------
  cat::MuonCollection selectedMuons;
  cat::MuonCollection vetoMuons;

  Handle<cat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  if( !doLooseLepton_ ){
    for( unsigned int i = 0; i < muons->size() ; i++ ){
      const cat::Muon & muon = muons->at(i);
      if     ( IsSelectMuon(muon) ) selectedMuons.push_back(muon);
      else if( IsVetoMuon(muon) )   vetoMuons.push_back(muon); // does not Include selected muons
    }
  }
  else if( !muons->empty() ){
    if( IsSelectMuon(muons->at(0)) ) selectedMuons.push_back( muons->at(0) );
    for( unsigned int i = 1; i < muons->size() ; i++ ){
      const cat::Muon & muon = muons->at(i);
      if( IsVetoMuon(muon) ) vetoMuons.push_back(muon);
    }
  }

  //---------------------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  TLorentzVector lepton;
  std::vector<float> *lepton_SF;
  lepton_SF = new std::vector<float>;

  int ch_tag  = 999;
  if( selectedMuons.size()     == 1 &&
      vetoMuons.size()         == 0 &&
      selectedElectrons.size() == 0 &&
      vetoElectrons.size()     == 0){
    lepton.SetPtEtaPhiE(selectedMuons[0].pt(), selectedMuons[0].eta(), selectedMuons[0].phi(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets
    b_Lepton_relIso = selectedMuons[0].relIso(0.4);
    b_Lepton_isIso = (b_Lepton_relIso < 0.15);

    if( isMC_ ){
      // Muon SF (Id/Iso/Trg)
      lepton_SF->push_back( SF_muonId_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()) ) );        // [0]-> IdSF
      lepton_SF->push_back( SF_muonId_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()),  1.0 ) );  // [1]-> IdSF+Error
      lepton_SF->push_back( SF_muonId_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()), -1.0 ) );  // [2]-> IdSF-Error
      lepton_SF->push_back( SF_muonIso_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()) ) );       // [3]-> IsoSF
      lepton_SF->push_back( SF_muonIso_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()),  1.0 ) ); // [4]-> IsoSF+Error
      lepton_SF->push_back( SF_muonIso_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()), -1.0 ) ); // [5]-> IsoSF-Error
      lepton_SF->push_back( SF_muonTrg_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()) ) );       // [6]-> TrgSF
      lepton_SF->push_back( SF_muonTrg_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()),  1.0 ) ); // [7]-> TrgSF+Error
      lepton_SF->push_back( SF_muonTrg_( selectedMuons[0].pt(), std::abs(selectedMuons[0].eta()), -1.0 ) ); // [8]-> TrgSF-Error
    }
    b_Electron_Scale->push_back( 1.0 ); // Dummy for muon
    b_Electron_Scale->push_back( 0.0 );
    b_Electron_Scale->push_back( 0.0 );
    b_Electron_Scale->push_back( 0.0 );
    b_Electron_Scale->push_back( 0.0 );
  }

  if( selectedMuons.size()     == 0 &&
      vetoMuons.size()         == 0 &&
      selectedElectrons.size() == 1 &&
      vetoElectrons.size()     == 0){
    ch_tag = 1; //electron + jets
    lepton.SetPtEtaPhiE(selectedElectrons[0].pt(), selectedElectrons[0].eta(), selectedElectrons[0].phi(), selectedElectrons[0].energy());
    b_Lepton_relIso = selectedElectrons.at(0).relIso(0.3);
    if( std::abs(selectedElectrons[0].scEta()) <= 1.479 and b_Lepton_relIso < (0.0287 + 0.506/static_cast<double>(lepton.Pt())) ) b_Lepton_isIso = true;
    else if( std::abs(selectedElectrons[0].scEta()) > 1.479 and b_Lepton_relIso < (0.0445 + 0.963/static_cast<double>(lepton.Pt())) ) b_Lepton_isIso = true;

    if( isMC_ ){
      // Lepton SF (Id/Reco/Zvtx)
      lepton_SF->push_back( SF_elecId_( selectedElectrons[0].pt(), selectedElectrons[0].scEta() ) );         // [0]-> IdSF
      lepton_SF->push_back( SF_elecId_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 ) );   // [1]-> IdSF+Error
      lepton_SF->push_back( SF_elecId_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 ) );   // [2]-> IdSF-Error
      lepton_SF->push_back( SF_elecReco_( selectedElectrons[0].pt(), selectedElectrons[0].scEta() ) );       // [3]-> RecoSF
      lepton_SF->push_back( SF_elecReco_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 ) ); // [4]-> RecoSF+Error
      lepton_SF->push_back( SF_elecReco_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 ) ); // [5]-> RecoSF-Error
      lepton_SF->push_back( SF_elecZvtx_( selectedElectrons[0].pt(), selectedElectrons[0].scEta() ) );       // [6]-> ZvtxSF
      lepton_SF->push_back( SF_elecZvtx_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 ) ); // [7]-> ZvtxSF+Error
      lepton_SF->push_back( SF_elecZvtx_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 ) ); // [8]-> ZvtxSF-Error
      lepton_SF->push_back( SF_elecTrg_( selectedElectrons[0].pt(), selectedElectrons[0].scEta() ) );        // [9]-> TrgSF, ttH v5
      lepton_SF->push_back( SF_elecTrg_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(),  1.0 ) );  // [10]-> TrgSF+Error
      lepton_SF->push_back( SF_elecTrg_( selectedElectrons[0].pt(), selectedElectrons[0].scEta(), -1.0 ) );  // [11]-> TrgSF-Error
      //b_Lepton_scEta = selectedElectrons[0].scEta();
    }
    b_Electron_Scale->push_back( selectedElectrons[0].smearedScale() );   // [0]-> Scale and smearing
    b_Electron_Scale->push_back( selectedElectrons[0].shiftedEnUp(0) );   // [1]-> energyScaleUp
    b_Electron_Scale->push_back( selectedElectrons[0].shiftedEnUp(1) );   // [2]-> energySigmaUp
    b_Electron_Scale->push_back( selectedElectrons[0].shiftedEnDown(0) ); // [3]-> energyScaleDown
    b_Electron_Scale->push_back( selectedElectrons[0].shiftedEnDown(1) ); // [4]-> energySigmaDown
  }

  //---------------------------------------------------------------------------
  // HLTrigger
  //---------------------------------------------------------------------------
  bool EvTrigger = false;
  bool IsTriggerMu = false;
  bool IsTriggerEl = false;
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
      (ch_tag == 1 && (IsTriggerEl || IsTriggerElHT) ) ){
    EvTrigger = true;
  }

  //---------------------------------------------------------------------------
  // Fill Tree with events that have ONLY one lepton
  //---------------------------------------------------------------------------
  // Check Gen Level for ttbar sample
  if( TTbarMC_ == 1 ){
    // ttbar Category
    if( !IsCat ) ch_tag = 999;
  }
  if( b_nGoodPV < 1 ) ch_tag = 999; 

  if( ch_tag<2 && EvTrigger && !METfiltered ){ // Single lepton event
    b_Channel  = ch_tag;

    b_Lepton_pt  = lepton.Pt();
    b_Lepton_eta = lepton.Eta();
    b_Lepton_phi = lepton.Phi();
    b_Lepton_e   = lepton.E();
    b_Lepton_SF  = lepton_SF;

    //---------------------------------------------------------------------------
    // Jets
    //---------------------------------------------------------------------------
    Handle<cat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    // Arrange Jets by CSV discriminator
    std::vector<int> JetIndex;
    for( unsigned int i = 0; i < jets->size() ; i++ ) JetIndex.push_back(i);

    int N_GoodJets = 0;
    int N_BJetsM = 0;

    // Initialize SF_btag
    // Jet_SF_CSV[Scenario][SystVariations];
    float Jet_SF_deepCSV[1][19];
    for( unsigned int ipTj=0; ipTj<1; ipTj++ ){
      for( unsigned int iu=0; iu<19; iu++ ){
        Jet_SF_deepCSV[ipTj][iu] = 1.0;
      }
    }

    int good_jet_count = 0;

    // Run again over all Jets
    for( unsigned int i = 0; i < JetIndex.size() ; i++ ){
      const cat::Jet & jet = jets->at(JetIndex[i]);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection (pT>20GeV to take into account SYST Variations)
      if( std::abs(jet.eta()) < 2.4 && jet.pt() > 20. && jet.tightLepVetoJetID() ){
        goodJet = true;
        good_jet_count += 1;
      }
      // Jet Cleaning
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if( dr_LepJet > 0.4 ) cleanJet = true;

      if( goodJet && cleanJet ){
        // Basic variables
        b_Jet_pt ->push_back(jet.pt());
        b_Jet_eta->push_back(jet.eta());
        b_Jet_phi->push_back(jet.phi());
        b_Jet_e  ->push_back(jet.energy());
        b_Jet_Index ->push_back(JetIndex[i]);
	
	      // Number of Jets (easy cross check)
        if( jet.pt() > 30. ) N_GoodJets++;

        // Parton Flavour
        b_Jet_partonFlavour->push_back(jet.partonFlavour());
        b_Jet_hadronFlavour->push_back(jet.hadronFlavour());

        // b-tag discriminant
        float jet_btagDis_deepCSVb = jet.bDiscriminator(BTAG_DeepCSVb);
        float jet_btagDis_deepCSVbb = jet.bDiscriminator(BTAG_DeepCSVbb);
        b_Jet_deepCSV ->push_back(jet_btagDis_deepCSVb+jet_btagDis_deepCSVbb);
        // c-tag discriminant
        float jet_btagDis_deepCvsL = jet.bDiscriminator(CTAG_DeepCSVc)/(jet.bDiscriminator(CTAG_DeepCSVc)+jet.bDiscriminator(CTAG_DeepCSVudsg));
        b_Jet_deepCvsL ->push_back(jet_btagDis_deepCvsL);
        float jet_btagDis_deepCvsB = jet.bDiscriminator(CTAG_DeepCSVc)/(jet.bDiscriminator(CTAG_DeepCSVc)+jet.bDiscriminator(BTAG_DeepCSVb)+jet.bDiscriminator(BTAG_DeepCSVbb));
        b_Jet_deepCvsB ->push_back(jet_btagDis_deepCvsB);
        if( jet_btagDis_deepCSVb+jet_btagDis_deepCSVbb > WP_BTAG_DeepCSVM ) N_BJetsM++;

        if( isMC_ ){
          b_Jet_JES_Up  ->push_back(jet.shiftedEnUp());
          b_Jet_JES_Down->push_back(jet.shiftedEnDown());

          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
          b_Jet_JER_Up   ->push_back( jer_valid(jet.smearedResUp()) );
          b_Jet_JER_Nom  ->push_back( jer_valid(jet.smearedRes()) );
          b_Jet_JER_Down ->push_back( jer_valid(jet.smearedResDown()) );

          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
          // Saving the central SF and the 18 syst. unc. for:
          if( jet.pt() > 30. ) for( unsigned int iu=0; iu<19; iu++ ) Jet_SF_deepCSV[0][iu] *= SF_deepCSV_.getSF(jet, iu);
        } // if(isMC_)	
      }// if(GoodJets)
    }// for(AllJets)
    
    // Number of Jets (easy cross check)
    b_Jet_NJet = N_GoodJets;
    b_Jet_NBJetM = N_BJetsM;

    for( unsigned int iu=0; iu<19; iu++ ) b_Jet_SF_deepCSV_30->push_back(1.0);
    (*b_Jet_SF_deepCSV_30)[0] = Jet_SF_deepCSV[0][0]; //Central
    // To save only the error
    for( unsigned int iu=1; iu<19; iu++ ){
      (*b_Jet_SF_deepCSV_30)[iu] = std::abs(Jet_SF_deepCSV[0][iu] - Jet_SF_deepCSV[0][0]); // Syst. Unc.
    }
   
    tree->Fill();
        
  }// if(ch_tag)
}

//------------- Good Muon Selection -----------------------
bool fcncLepJetsAnalyzer::IsSelectMuon(const cat::Muon & i_muon_candidate){

  bool GoodMuon=true;

  // Tight selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoTight));
  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 20);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  return GoodMuon;
}

//------------- Loose Muon Selection -----------------------
bool fcncLepJetsAnalyzer::IsVetoMuon(const cat::Muon & i_muon_candidate){

  bool GoodMuon=true;

  // Loose selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.passed(reco::Muon::CutBasedIdLoose|reco::Muon::PFIsoLoose));
  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 15);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool fcncLepJetsAnalyzer::IsSelectElectron(const cat::Electron & i_electron_candidate){

  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // Electron cut based selection, wp80 is tight
  // https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_Recipe_for_regul
  if( !doLooseLepton_ ) GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-tight") > 0.0;
  else                  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-tight-noiso") > 0.0;

  return GoodElectron;

}

//------------- Loose Electron Selection -----------------------
bool fcncLepJetsAnalyzer::IsVetoElectron(const cat::Electron & i_electron_candidate){

  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_Recipe_for_regul
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-veto") > 0.0;

  return GoodElectron;
}

float fcncLepJetsAnalyzer::weight_valid(float weight){
  if( std::isnan(weight) or std::isinf(weight) ) return 1.0;
  if( (weight < 0.01) or (weight > 100) )        return 1.0;
  return weight;
};

float fcncLepJetsAnalyzer::jer_valid(float weight){
  if( std::isnan(weight) or std::isinf(weight) ) return 0.0;
  if( weight < 0.0 )                             return 0.0;
  return weight;
};

//define this as a plug-in
DEFINE_FWK_MODULE(fcncLepJetsAnalyzer);
