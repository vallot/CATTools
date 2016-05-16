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

  bool isMC_;
  // TTbarMC_ == 0, No ttbar
  // TTbarMC_ == 1, ttbar Signal
  // TTbarMC_ == 2, ttbar Background
  int TTbarMC_;

  edm::EDGetTokenT<float>                        genWeightToken_;
  edm::EDGetTokenT<vector<float>>                pdfWeightsToken_;
  edm::EDGetTokenT<vector<float>>                scaleUpWeightsToken_;
  edm::EDGetTokenT<vector<float>>                scaleDownWeightsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>  genToken_;
  edm::EDGetTokenT<reco::GenJetCollection>       genJetToken_;
  edm::EDGetTokenT<cat::GenTopCollection>        genttbarCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
  edm::EDGetTokenT<cat::JetCollection>           jetToken_;
  edm::EDGetTokenT<cat::METCollection>           metToken_;
  edm::EDGetTokenT<int>                          pvToken_;
  edm::EDGetTokenT<float>                        puWeightToken_;
  edm::EDGetTokenT<float>                        puUpWeightToken_;
  edm::EDGetTokenT<float>                        puDownWeightToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
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
  std::vector<int> *b_GenConeCatID;

  // MET
  float b_MET, b_MET_phi;

  // GEN Leptons
  float b_GenLepton_pT;
  float b_GenLepton_eta;
  // Leptons
  float b_Lepton_px;
  float b_Lepton_py;
  float b_Lepton_pz;
  float b_Lepton_E;
  float b_Lepton_pT;
  std::vector<float> *b_Lepton_SF;
  float b_Lepton_LES;

  // GEN Jets
  std::vector<float> *b_GenJet_pT;
  // Jets
  int b_Jet_Number;
  std::vector<float> *b_Jet_px;
  std::vector<float> *b_Jet_py;
  std::vector<float> *b_Jet_pz;
  std::vector<float> *b_Jet_E;
  std::vector<float> *b_Jet_pT;
  // Flavour
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
  std::vector<float> *b_Jet_SF_CSV;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Histograms
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Number of events
  TH1F *EventInfo;

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
  isMC_    (iConfig.getUntrackedParameter<bool>("sampleLabel",  true)),
  TTbarMC_ (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0))
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

  //now do what ever initialization is needed
  genWeightToken_       = consumes<float>(iConfig.getParameter<edm::InputTag>("genWeightLabel"));
  pdfWeightsToken_       = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("pdfWeightLabel"));
  scaleUpWeightsToken_   = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaleUpWeightLabel"));
  scaleDownWeightsToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaleDownWeightLabel"));

  genToken_             = consumes<reco::GenParticleCollection>  (iConfig.getParameter<edm::InputTag>("genLabel"));
  genJetToken_          = consumes<reco::GenJetCollection>       (iConfig.getParameter<edm::InputTag>("genJetLabel"));
  genttbarCatToken_     = consumes<cat::GenTopCollection>        (iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));

  muonToken_         = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_     = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_          = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_          = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_           = consumes<int>                          (iConfig.getParameter<edm::InputTag>("pvLabel"));

  puWeightToken_     = consumes<float>                        (iConfig.getParameter<edm::InputTag>("puWeightLabel"));
  puUpWeightToken_   = consumes<float>                        (iConfig.getParameter<edm::InputTag>("puUpWeightLabel"));
  puDownWeightToken_ = consumes<float>                        (iConfig.getParameter<edm::InputTag>("puDownWeightLabel"));

  triggerBits_       = consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_    = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerObjects"));


  b_PUWeight    = new std::vector<float>;
  b_ScaleWeight = new std::vector<float>;

  b_GenConeCatID = new std::vector<int>;

  b_GenJet_pT = new std::vector<float>;

  b_Lepton_SF  = new std::vector<float>;

  b_Jet_px   = new std::vector<float>;
  b_Jet_py   = new std::vector<float>;
  b_Jet_pz   = new std::vector<float>;
  b_Jet_E    = new std::vector<float>;
  b_Jet_pT   = new std::vector<float>;
  b_Jet_CSV  = new std::vector<float>;

  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;

  b_Jet_JES_Up   = new std::vector<float>;
  b_Jet_JES_Down = new std::vector<float>;
  b_Jet_JER_Up   = new std::vector<float>;
  b_Jet_JER_Nom  = new std::vector<float>;
  b_Jet_JER_Down = new std::vector<float>;

  b_Jet_SF_CSV  = new std::vector<float>;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "Event/I");
  tree->Branch("run",        &b_Run,         "Run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "Lumi_Number/I");
  tree->Branch("genweight",  &b_GenWeight,   "GenWeight/F");

  tree->Branch("PUWeight", "std::vector<float>", &b_PUWeight);
  tree->Branch("GoodPV",   &b_nGoodPV,  "nGoodPV/I");

  tree->Branch("channel",  &b_Channel,  "Channel/I");

  tree->Branch("MET",     &b_MET,     "MET/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("lepton_px", &b_Lepton_px, "lepton_px/F");
  tree->Branch("lepton_py", &b_Lepton_py, "lepton_py/F");
  tree->Branch("lepton_pz", &b_Lepton_pz, "lepton_pz/F");
  tree->Branch("lepton_E" , &b_Lepton_E,  "lepton_E/F" );
  tree->Branch("lepton_pT", &b_Lepton_pT, "lepton_pT/F" );

  tree->Branch("lepton_SF",  "std::vector<float>", &b_Lepton_SF );
  tree->Branch("lepton_LES", &b_Lepton_LES, "lepton_LES/F" );

  tree->Branch("jet_px", "std::vector<float>", &b_Jet_px);
  tree->Branch("jet_py", "std::vector<float>", &b_Jet_py);
  tree->Branch("jet_pz", "std::vector<float>", &b_Jet_pz);
  tree->Branch("jet_E" , "std::vector<float>", &b_Jet_E );
  tree->Branch("jet_pT", "std::vector<float>", &b_Jet_pT );

  tree->Branch("jet_CSV" ,   "std::vector<float>", &b_Jet_CSV );
  tree->Branch("jet_SF_CSV", "std::vector<float>", &b_Jet_SF_CSV );

  tree->Branch("jet_Number" , &b_Jet_Number, "jet_number/I" );

  tree->Branch("jet_partonFlavour", "std::vector<int>", &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>", &b_Jet_hadronFlavour);

  tree->Branch("jet_JES_Up",   "std::vector<float>", &b_Jet_JES_Up );
  tree->Branch("jet_JES_Down", "std::vector<float>", &b_Jet_JES_Down );
  tree->Branch("jet_JER_Up",   "std::vector<float>", &b_Jet_JER_Up );
  tree->Branch("jet_JER_Nom",  "std::vector<float>", &b_Jet_JER_Nom );
  tree->Branch("jet_JER_Down", "std::vector<float>", &b_Jet_JER_Down );


  // GEN Tree (only ttbarSignal)
  if(isMC_ && TTbarMC_==1){
    tree->Branch("scaleweight", "std::vector<float>", &b_ScaleWeight );
    tree->Branch("genconecatid" , "std::vector<int>", &b_GenConeCatID);
    tree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    tree->Branch("genlepton_pT",  &b_GenLepton_pT,  "genlepton_pT/F");
    tree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    tree->Branch("genjet_pT", "std::vector<float>", &b_GenJet_pT);

    gentree = fs->make<TTree>("gentree", "TopGENTree");
    gentree->Branch("genchannel",   &b_GenChannel,   "genchannel/I");
    gentree->Branch("genconecatid" , "std::vector<int>", &b_GenConeCatID);
    gentree->Branch("genlepton_pT", &b_GenLepton_pT, "genlepton_pT/F");
    gentree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    gentree->Branch("genjet_pT", "std::vector<float>", &b_GenJet_pT);
  }

  EventInfo = fs->make<TH1F>("EventInfo","Event Information",14,0,14);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
  EventInfo->GetXaxis()->SetBinLabel(3,"ttbb Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(4,"ttbb Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(5,"ttbb Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(6,"ttjj Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(7,"ttjj Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(8,"ttjj Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(9, "muR=Nom  muF=Up");
  EventInfo->GetXaxis()->SetBinLabel(10,"muR=Nom  muF=Down");
  EventInfo->GetXaxis()->SetBinLabel(11,"muR=Up  muF=Nom");
  EventInfo->GetXaxis()->SetBinLabel(12,"muR=Up  muF=Up");
  EventInfo->GetXaxis()->SetBinLabel(13,"muR=Down  muF=Nom");
  EventInfo->GetXaxis()->SetBinLabel(14,"muR=Down  muF=Down");
}


ttbbLepJetsAnalyzer::~ttbbLepJetsAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete b_PUWeight;
  delete b_ScaleWeight;

  delete b_GenConeCatID;

  delete b_GenJet_pT;

  delete b_Lepton_SF;

  delete b_Jet_px;
  delete b_Jet_py;
  delete b_Jet_pz;
  delete b_Jet_E;
  delete b_Jet_pT;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_JES_Up;
  delete b_Jet_JES_Down;
  delete b_Jet_JER_Up;
  delete b_Jet_JER_Nom;
  delete b_Jet_JER_Down;

  delete b_Jet_SF_CSV;

  delete b_Jet_CSV;

}

//
// member functions
//

// ------------ method called for each event  ------------
void ttbbLepJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_PUWeight->clear();
  b_ScaleWeight->clear();

  b_GenConeCatID->clear();

  b_GenJet_pT->clear();

  b_Lepton_SF->clear();

  b_Jet_px->clear();
  b_Jet_py->clear();
  b_Jet_pz->clear();
  b_Jet_E->clear();
  b_Jet_pT->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_JES_Up->clear();
  b_Jet_JES_Down->clear();
  b_Jet_JER_Up->clear();
  b_Jet_JER_Nom->clear();
  b_Jet_JER_Down->clear();

  b_Jet_CSV->clear();
  b_Jet_SF_CSV->clear();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

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
    edm::Handle<vector<float>> scaleUpWeightsHandle, scaleDownWeightsHandle;
    iEvent.getByToken(scaleUpWeightsToken_, scaleUpWeightsHandle);
    iEvent.getByToken(scaleDownWeightsToken_, scaleDownWeightsHandle);

    // muR/muF Scale Weights
    // Comment from J. Goh - scale up/down are split in the new format, but I'm keeping the original index to minimize change
    // Originally, it was (Nom,Up) (Nom,Down), (Up,Nom), (Up, Up), (Down,Nom), (Down,Down)
    b_ScaleWeight->push_back(scaleUpWeightsHandle->at(0));
    b_ScaleWeight->push_back(scaleDownWeightsHandle->at(0));
    b_ScaleWeight->push_back(scaleUpWeightsHandle->at(1));
    b_ScaleWeight->push_back(scaleUpWeightsHandle->at(2));
    b_ScaleWeight->push_back(scaleDownWeightsHandle->at(1));
    b_ScaleWeight->push_back(scaleDownWeightsHandle->at(2));

    // Sum of muR/muF Scale Weights
    EventInfo->Fill(8.5,  b_ScaleWeight->at(0));
    EventInfo->Fill(9.5,  b_ScaleWeight->at(1));
    EventInfo->Fill(10.5, b_ScaleWeight->at(2));
    EventInfo->Fill(11.5, b_ScaleWeight->at(3));
    EventInfo->Fill(12.5, b_ScaleWeight->at(4));
    EventInfo->Fill(13.5, b_ScaleWeight->at(5));
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Generated Particles (For Pythia8)
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  int nGenLep  = 999;

  if(TTbarMC_ > 0) {

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

    if(isMC_ && TTbarMC_==1){
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

	edm::Handle<reco::GenJetCollection> genJets;
	iEvent.getByToken(genJetToken_, genJets);
	
	for (unsigned int j = 0; j < genJets->size(); j++){
	  const cat::GenJet & gjet = (*genJets)[j];
	  if(std::abs(gjet.eta())< 2.5 &&
	     gjet.pt()> 20){
	    b_GenJet_pT->push_back(gjet.pt());
	  }// if(Good genJet)
	}// for(genJet)	
	
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
    lepton.SetPxPyPzE(selectedMuons[0].px(), selectedMuons[0].py(), selectedMuons[0].pz(), selectedMuons[0].energy());
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
    lepton.SetPxPyPzE(selectedElectrons[0].px(), selectedElectrons[0].py(), selectedElectrons[0].pz(), selectedElectrons[0].energy());

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

  bool EvTrigger = false;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  if ( (ch_tag == 0 && (trigHelper.triggerFired("HLT_IsoMu20_v") || trigHelper.triggerFired("HLT_IsoTkMu20_v"))) ||
       (ch_tag == 1 && (trigHelper.triggerFired("HLT_Ele22_eta2p1_WP75_Gsf_v") || trigHelper.triggerFired("HLT_Ele22_eta2p1_WPLoose_Gsf_v"))) ) {
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
  } // if(TTbarMC_ >0)

  if (ch_tag<2 && EvTrigger){ // Single lepton event

    b_Channel  = ch_tag;

    b_Lepton_px = lepton.Px();
    b_Lepton_py = lepton.Py();
    b_Lepton_pz = lepton.Pz();
    b_Lepton_E  = lepton.E();
    b_Lepton_pT = lepton.Pt();

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
    float Jet_SF_CSV[19];
    for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] = 1.0;

    // Run again over all Jets (CSV order)
    for (unsigned int i = 0; i < JetIndex.size() ; i++) {

      const cat::Jet & jet = jets->at(JetIndex[i]);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection (pT>15GeV to take into account SYST Variations)
      if(std::abs(jet.eta()) < 2.4 && jet.pt() > 15. && jet.LooseId()) goodJet = true;
      // Jet Cleaning
      TLorentzVector vjet(jet.px(), jet.py(), jet.pz(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if(dr_LepJet > 0.4) cleanJet = true;

      if(goodJet && cleanJet){
        // Basic variables
        N_GoodJets ++;
        b_Jet_px->push_back(jet.px());
        b_Jet_py->push_back(jet.py());
        b_Jet_pz->push_back(jet.pz());
        b_Jet_E ->push_back(jet.energy());
        b_Jet_pT->push_back(jet.pt());

        // Parton Flavour
        b_Jet_partonFlavour->push_back(jet.partonFlavour());
        b_Jet_hadronFlavour->push_back(jet.hadronFlavour());

        // b-tag discriminant
        float jet_btagDis_CSV = jet.bDiscriminator(BTAG_CSVv2);
        b_Jet_CSV ->push_back(jet_btagDis_CSV);

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
          if(jet.pt() > 25.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] *= SF_CSV_.getSF(jet, iu);
        } // if(isMC_)
	
      }
    }

    b_Jet_Number = N_GoodJets;

    for (unsigned int iu=0; iu<19; iu++) b_Jet_SF_CSV->push_back(1.0);
    (*b_Jet_SF_CSV)[0] = Jet_SF_CSV[0]; //Central
    // To save only the error
    for (unsigned int iu=1; iu<19; iu+=2) (*b_Jet_SF_CSV)[iu] = Jet_SF_CSV[iu] - Jet_SF_CSV[0] ; // Syst. Unc. Up
    for (unsigned int iu=2; iu<19; iu+=2) (*b_Jet_SF_CSV)[iu] = Jet_SF_CSV[0]  - Jet_SF_CSV[iu]; // Syst. Unc. Down

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Fill Tree with event at 1 lepton cut level
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    tree->Fill();

  } // if(ch_tag)

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
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0.0;

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
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto") > 0;

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

