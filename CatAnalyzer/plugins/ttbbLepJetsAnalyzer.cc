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

#include "CATTools/CommonTools/interface/AnalysisHelper.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// b-tagging Eff
//#include "CATTools/CatAnalyzer/interface/BTagSFUtil.h"

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
  bool IsTightMuon    (const cat::Muon     & i_muon_candidate);
  bool IsLooseMuon    (const cat::Muon     & i_muon_candidate);
  bool IsTightElectron(const cat::Electron & i_electron_candidate);
  bool IsLooseElectron(const cat::Electron & i_electron_candidate);

  bool isMC_;
  // TTbarMC_ == 0, No ttbar
  // TTbarMC_ == 1, ttbar Signal
  // TTbarMC_ == 2, ttbar Background
  int TTbarMC_;

  edm::EDGetTokenT<float>                        genWeightToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>  genToken_;
  edm::EDGetTokenT<reco::GenJetCollection>       genJetToken_;
  edm::EDGetTokenT<int>                          genttbarCatToken_;
  edm::EDGetTokenT<cat::GenTopCollection>        genttbarConeCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
  edm::EDGetTokenT<cat::JetCollection>           jetToken_;
  edm::EDGetTokenT<cat::METCollection>           metToken_;
  edm::EDGetTokenT<int>                          pvToken_;
  edm::EDGetTokenT<float>                        puWeightToken_;

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

  // PU/Vertices
  float b_PUWeight;
  int b_nGoodPV;

  // Channel and Categorization
  int b_Channel;
  int b_GenCatID;
  std::vector<int> *b_GenConeCatID;
  // GEN Channel
  int b_GenChannel;

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

  // GEN Jets
  std::vector<float> *b_GenJet_pT;
  std::vector<float> *b_GenJet_Flavour;
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
  // Smearing and Shifts
  std::vector<float> *b_Jet_smearedRes;
  std::vector<float> *b_Jet_smearedResDown;
  std::vector<float> *b_Jet_smearedResUp;
  std::vector<float> *b_Jet_shiftedEnUp;
  std::vector<float> *b_Jet_shiftedEnDown;
  // b-Jet discriminant
  std::vector<float> *b_Jet_CSV;
  std::vector<float> *b_Jet_SVM;


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Histograms
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Number of events
  TH1F *EventInfo;

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
  //now do what ever initialization is needed
  genWeightToken_       = consumes<float>                        (iConfig.getParameter<edm::InputTag>("genWeightLabel"));
  genToken_             = consumes<reco::GenParticleCollection>  (iConfig.getParameter<edm::InputTag>("genLabel"));
  genJetToken_          = consumes<reco::GenJetCollection>       (iConfig.getParameter<edm::InputTag>("genJetLabel"));
  genttbarCatToken_     = consumes<int>                          (iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));
  genttbarConeCatToken_ = consumes<cat::GenTopCollection>        (iConfig.getParameter<edm::InputTag>("genttbarConeCatLabel"));

  muonToken_         = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_     = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_          = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_          = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_           = consumes<int>                          (iConfig.getParameter<edm::InputTag>("pvLabel"));
  puWeightToken_     = consumes<float>                        (iConfig.getParameter<edm::InputTag>("puWeightLabel"));

  triggerBits_       = consumes<edm::TriggerResults>          (iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_    = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerObjects"));

  b_GenConeCatID= new std::vector<int>;

  b_GenJet_pT      = new std::vector<float>;
  b_GenJet_Flavour = new std::vector<float>;

  b_Jet_px = new std::vector<float>;
  b_Jet_py = new std::vector<float>;
  b_Jet_pz = new std::vector<float>;
  b_Jet_E  = new std::vector<float>;
  b_Jet_pT = new std::vector<float>;

  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;

  b_Jet_smearedRes     = new std::vector<float>;
  b_Jet_smearedResDown = new std::vector<float>;
  b_Jet_smearedResUp   = new std::vector<float>;
  b_Jet_shiftedEnUp    = new std::vector<float>;
  b_Jet_shiftedEnDown  = new std::vector<float>;

  b_Jet_CSV  = new std::vector<float>;
  b_Jet_SVM  = new std::vector<float>;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "Event/I");
  tree->Branch("run",        &b_Run,         "Run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "Lumi_Number/I");
  tree->Branch("genweight",  &b_GenWeight,   "GenWeight/F");

  tree->Branch("PUWeight", &b_PUWeight, "PUWeight/F");
  tree->Branch("GoodPV",   &b_nGoodPV,  "nGoodPV/I");

  tree->Branch("channel",  &b_Channel,  "Channel/I");

  tree->Branch("MET",     &b_MET,     "MET/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("lepton_px", &b_Lepton_px, "lepton_px/F");
  tree->Branch("lepton_py", &b_Lepton_py, "lepton_py/F");
  tree->Branch("lepton_pz", &b_Lepton_pz, "lepton_pz/F");
  tree->Branch("lepton_E" , &b_Lepton_E,  "lepton_E/F" );
  tree->Branch("lepton_pT", &b_Lepton_pT, "lepton_pT/F" );

  tree->Branch("jet_px", "std::vector<float>", &b_Jet_px);
  tree->Branch("jet_py", "std::vector<float>", &b_Jet_py);
  tree->Branch("jet_pz", "std::vector<float>", &b_Jet_pz);
  tree->Branch("jet_E" , "std::vector<float>", &b_Jet_E );
  tree->Branch("jet_pT", "std::vector<float>", &b_Jet_pT );

  tree->Branch("jet_Number" , &b_Jet_Number, "jet_number/I" );

  tree->Branch("jet_partonFlavour", "std::vector<int>", &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>", &b_Jet_hadronFlavour);

  tree->Branch("jet_smearedRes",     "std::vector<float>", &b_Jet_smearedRes);
  tree->Branch("jet_smearedResDown", "std::vector<float>", &b_Jet_smearedResDown);
  tree->Branch("jet_smearedResUp",   "std::vector<float>", &b_Jet_smearedResUp);
  tree->Branch("jet_shiftedEnUp",    "std::vector<float>", &b_Jet_shiftedEnUp);
  tree->Branch("jet_shiftedEnDown",  "std::vector<float>", &b_Jet_shiftedEnDown);

  tree->Branch("jet_CSV" , "std::vector<float>", &b_Jet_CSV );
  tree->Branch("jet_SVM" , "std::vector<float>", &b_Jet_SVM );

  // GEN Tree (only ttbarSignal)
  if(isMC_ && TTbarMC_==1){
    tree->Branch("gencatid", &b_GenCatID, "GenCatID/I");
    tree->Branch("genconecatid" , "std::vector<int>", &b_GenConeCatID);
    tree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    tree->Branch("genlepton_pT",  &b_GenLepton_pT,  "genlepton_pT/F");
    tree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    tree->Branch("genjet_pT", "std::vector<float>", &b_GenJet_pT);
    tree->Branch("genjet_Flavour", "std::vector<float>", &b_GenJet_Flavour);
    
    gentree = fs->make<TTree>("gentree", "TopGENTree");
    gentree->Branch("genchannel",   &b_GenChannel,   "genchannel/I");
    gentree->Branch("gencatid", &b_GenCatID, "GenCatID/I");
    gentree->Branch("genconecatid" , "std::vector<int>", &b_GenConeCatID);
    gentree->Branch("genlepton_pT", &b_GenLepton_pT, "genlepton_pT/F");
    gentree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    gentree->Branch("genjet_pT", "std::vector<float>", &b_GenJet_pT);
    gentree->Branch("genjet_Flavour", "std::vector<float>", &b_GenJet_Flavour);
  }

  EventInfo = fs->make<TH1F>("EventInfo","Event Information",5,0,5); 
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
  EventInfo->GetXaxis()->SetBinLabel(3,"Number of ttbb Events");
  EventInfo->GetXaxis()->SetBinLabel(4,"Sum of ttbb Weights");

}


ttbbLepJetsAnalyzer::~ttbbLepJetsAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete b_GenConeCatID;

  delete b_GenJet_pT;
  delete b_GenJet_Flavour;

  delete b_Jet_px;
  delete b_Jet_py;
  delete b_Jet_pz;
  delete b_Jet_E;
  delete b_Jet_pT;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_smearedRes;
  delete b_Jet_smearedResDown;
  delete b_Jet_smearedResUp;
  delete b_Jet_shiftedEnUp;
  delete b_Jet_shiftedEnDown;

  delete b_Jet_CSV;
  delete b_Jet_SVM;

}

//
// member functions
//

// ------------ method called for each event  ------------
void ttbbLepJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_GenConeCatID->clear();
  b_GenJet_pT->clear();
  b_GenJet_Flavour->clear();

  b_Jet_px->clear();
  b_Jet_py->clear();
  b_Jet_pz->clear();
  b_Jet_E->clear();
  b_Jet_pT->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_smearedRes->clear();
  b_Jet_smearedResDown->clear();
  b_Jet_smearedResUp->clear();
  b_Jet_shiftedEnUp->clear();
  b_Jet_shiftedEnDown->clear();

  b_Jet_CSV->clear();
  b_Jet_SVM->clear();

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
    b_PUWeight = *PUWeight;

    //---------------------------------------------------------------------------
    // Weights at Generation Level: MC@NLO
    //---------------------------------------------------------------------------
    
    edm::Handle<float> genWeight;
    
    iEvent.getByToken(genWeightToken_, genWeight);
    b_GenWeight = *genWeight;

    EventInfo->Fill(0.5, 1.0);         // First bin: Number of Events
    EventInfo->Fill(1.5, b_GenWeight); // Second bin: Sum of Weights
  }
  
  else{ 
    b_PUWeight = 1;
    b_GenWeight = 1.0;
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Generated Particles (For Pythia8)
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  int nGenLep  = 999;
  bool GenMu_p = false;
  bool GenMu_m = false;
  bool GenEG_p = false;
  bool GenEG_m = false;

  bool GenTau_p = false;
  bool GenTau_m = false;

  b_GenCatID = 0; // Preset for non ttbar samples

  if(TTbarMC_ > 0) {

    //---------------------------------------------------------------------------
    // Event Categorization
    // Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GenHFHadronMatcher
    //---------------------------------------------------------------------------
    edm::Handle<int> genttbarCatHandle;
    iEvent.getByToken( genttbarCatToken_, genttbarCatHandle );
    b_GenCatID = *genttbarCatHandle;

    int tu = *genttbarCatHandle % 100;

    if(tu == 55 || tu == 54 || tu == 53){
      EventInfo->Fill(2.5, 1.0);         // Third bin: Number of ttbb Events
      EventInfo->Fill(3.5, b_GenWeight); // Fourth bin: Sum of ttbb Weights
    }


    //---------------------------------------------------------------------------
    // Event Categorization Using Cone
    //---------------------------------------------------------------------------
    edm::Handle<cat::GenTopCollection> genttbarConeCat;
    iEvent.getByToken( genttbarConeCatToken_, genttbarConeCat );

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

    //---------------------------------------------------------------------------
    // GEN Particles: Decay mode
    //---------------------------------------------------------------------------

    TLorentzVector Genlepton;

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genToken_, genParticles);

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
      if(abs(id) == 11 || abs(id) == 13 || abs(id) == 15){

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
	  if (abs(id) == 13 || abs(id) == 11){
	    Genlepton.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy());  
	    if (id == 13)  GenMu_p = true; // Muons +
	    if (id == 11)  GenEG_p = true; // Electrons +
	    if (id == -13) GenMu_m = true; // Muons -
	    if (id == -11) GenEG_m = true; // Electrons -
	  }
	  if (abs(id) == 15){         // Taus 
            if (id == 15)   GenTau_p = true; // Muons +
            if (id == -15)  GenTau_m = true; // Electrons +
            for(unsigned int h = 0; h <  gp.numberOfDaughters(); h++) {
              const reco::Candidate *gd = gp.daughter(h);
              const int taudauid = gd->pdgId();
	      if (abs(taudauid) == 13 || abs(taudauid) == 11){
		Genlepton.SetPxPyPzE(gd->px(), gd->py(), gd->pz(), gd->energy());  
		if (taudauid == 13)  GenMu_p = true; // Muons +
		if (taudauid == 11)  GenEG_p = true; // Electrons +
		if (taudauid == -13) GenMu_m = true; // Muons -
		if (taudauid == -11) GenEG_m = true; // Electrons -
	      }
	    } // for(taus' daughters)
	  } // if(taus)
	  
        } // if(t->W)
	
      }// if (mu || e || tau)
    } // for(genParticles)

    // nGenLep = 0; Full Hadronic
    // nGenLep = 1; Single Lepton
    // nGenLep = 2; Dilepton
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !! Problems
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // if(!GenMu_p && !GenEG_p && !GenMu_m && !GenEG_m && !GenTau_m && !GenTau_p) nGenLep = 0;
    // else if((GenMu_p  && GenMu_m)  || 
    // 	    (GenEG_p  && GenEG_m)  ||
    // 	    (GenTau_p && GenTau_m) ||
    // 	    (GenMu_p  && GenEG_m)  ||
    // 	    (GenMu_p  && GenTau_m) ||
    // 	    (GenEG_p  && GenTau_m) ||
    // 	    (GenEG_p  && GenMu_m)  ||
    // 	    (GenTau_p && GenMu_m)  ||
    // 	    (GenTau_p && GenEG_m) ) nGenLep = 2;
    // else nGenLep = 1;

    // Using the GenChannel from GenTop categorization
    nGenLep = genttbarConeCat->begin()->semiLeptonic(0);

    if(!GenMu_p && !GenEG_p && !GenMu_m && !GenEG_m && !GenTau_m && !GenTau_p);

    if(isMC_ && TTbarMC_==1){
      if(nGenLep == 1){

	// if (GenMu_m || GenMu_p)      b_GenChannel = 0;
	// else if (GenEG_m || GenEG_p) b_GenChannel = 1;
	
	// b_GenLepton_pT  = Genlepton.Pt();
	// b_GenLepton_eta = Genlepton.Eta();
	
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
	  //const reco::Candidate & gjet = (*genJets)[j];	  
	  const cat::GenJet & gjet = (*genJets)[j];
	  if(std::abs(gjet.eta())< 2.5 && 
	     gjet.pt()> 20){
	    b_GenJet_pT->push_back(gjet.pt());
	    b_GenJet_Flavour->push_back(gjet.partonFlavour());
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

    if( IsTightElectron( electron ) ) selectedElectrons.push_back( electron );
    else if( IsLooseElectron( electron ) ) vetoElectrons.push_back( electron ); // does not Include selected electrons

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

    if( IsTightMuon( muon) ) selectedMuons.push_back( muon);
    else if( IsLooseMuon( muon) ) vetoMuons.push_back( muon); // does not Include selected muons

  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  TLorentzVector lepton;
  int ch_tag  = 999;

  if(selectedMuons.size()     == 1 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 0 &&
     vetoElectrons.size()     == 0){
    lepton.SetPxPyPzE(selectedMuons[0].px(), selectedMuons[0].py(), selectedMuons[0].pz(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets
  }

  if(selectedMuons.size()     == 0 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 1 &&
     vetoElectrons.size()     == 0){
    lepton.SetPxPyPzE(selectedElectrons[0].px(), selectedElectrons[0].py(), selectedElectrons[0].pz(), selectedElectrons[0].energy());
    ch_tag = 1; //electron + jets
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // HLTrigger
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  bool EvTrigger = false; // Trigger requirement not yet applied
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

        float iJetCSV = jet_i.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        float jJetCSV = jet_j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        if(jJetCSV > iJetCSV){
          float tempIndex = JetIndex[ijet];
          JetIndex[ijet] = JetIndex[jjet];
          JetIndex[jjet] = tempIndex;
        }//if(jJetCSV > iJetCSV)
      }
    }

    int N_GoodJets = 0;
    for (unsigned int i = 0; i < JetIndex.size() ; i++) {
      
      const cat::Jet & jet = jets->at(JetIndex[i]);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection
      if(std::abs(jet.eta()) < 2.4 && jet.pt() > 25. && jet.LooseId()) goodJet = true;
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

        // Smeared and Shifted
        b_Jet_smearedRes     ->push_back(jet.smearedRes() );
        b_Jet_smearedResDown ->push_back(jet.smearedResDown());
        b_Jet_smearedResUp   ->push_back(jet.smearedResUp());
        b_Jet_shiftedEnUp    ->push_back(jet.shiftedEnUp());
        b_Jet_shiftedEnDown  ->push_back(jet.shiftedEnDown());

        b_Jet_SVM ->push_back(jet.vtxMass());

        // b-tag discriminant
        float jet_btagDis_CSV = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        b_Jet_CSV ->push_back(jet_btagDis_CSV);

      }
    }

    b_Jet_Number = N_GoodJets;
    // Fill Tree with event at 1 lepton cut level
    tree->Fill();

  } // if(ch_tag)

}

//------------- Good Muon Selection -----------------------
bool ttbbLepJetsAnalyzer::IsTightMuon(const cat::Muon & i_muon_candidate)
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
bool ttbbLepJetsAnalyzer::IsLooseMuon(const cat::Muon & i_muon_candidate)
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
bool ttbbLepJetsAnalyzer::IsTightElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.1);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // Electron cut based selection
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  // GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") > 0.0;

  // Electron MVA selection (Tight: WP80)
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_for_7_4_12_Spring15_MVA
  if (std::abs(i_electron_candidate.eta()) < 0.8)        GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.988153; // Barrel
  else if (std::abs(i_electron_candidate.eta()) > 0.8 && 
  	   std::abs(i_electron_candidate.eta()) < 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.967910; // Barrel
  else if (std::abs(i_electron_candidate.eta()) > 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.841729; // EndCaps

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;

  GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.0766 );

  // Effective Area Parametrization can be found in:
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}
//------------- Loose Electron Selection -----------------------
bool ttbbLepJetsAnalyzer::IsLooseElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto") > 0;

  // Electron MVA selection (Medium: WP90)
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_for_7_4_12_Spring15_MVA
  // if(std::abs(i_electron_candidate.eta()) < 0.8)        GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") > 0.972153;
  // else if(std::abs(i_electron_candidate.eta()) > 0.8 &&
  // 	  std::abs(i_electron_candidate.eta()) < 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") > 0.922126;
  // else if(std::abs(i_electron_candidate.eta()) > 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") > 0.610764;

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

