// -*- C++ -*-
//
// Package:    TtbarSingleLeptonAnalyzer
// Class:      TtbarSingleLeptonAnalyzer
//
/**\class TtbarSingleLeptonAnalyzer TtbarSingleLeptonAnalyzer.cc

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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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
#include "TFile.h"
#include "TLorentzVector.h"
//
// class declaration
//

using namespace cat;

class TtbarSingleLeptonAnalyzer : public edm::EDAnalyzer {
public:
  explicit TtbarSingleLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarSingleLeptonAnalyzer();

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

  edm::EDGetTokenT<edm::View<reco::GenParticle> >  genToken_;
  edm::EDGetTokenT<edm::View<cat::Muon> >          muonToken_;
  edm::EDGetTokenT<edm::View<cat::Electron> >      electronToken_;
  edm::EDGetTokenT<edm::View<cat::Jet> >           jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> >           metToken_;
  edm::EDGetTokenT<int>                            pvToken_;
  edm::EDGetTokenT<float>                          puWeight_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  // ----------member data ---------------------------

  TTree *tree;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Tree Branches
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Event info
  int b_Event, b_Run, b_Lumi_Number;

  // PU/Vertices
  float b_PUWeight;
  int b_nGoodPV;

  int b_Channel;

  // MET
  float b_MET, b_MET_phi;

  // Leptons
  float b_Lepton_px;
  float b_Lepton_py;
  float b_Lepton_pz;
  float b_Lepton_E;

  // Jets
  std::vector<float> *b_Jet_px;
  std::vector<float> *b_Jet_py;
  std::vector<float> *b_Jet_pz;
  std::vector<float> *b_Jet_E;
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
TtbarSingleLeptonAnalyzer::TtbarSingleLeptonAnalyzer(const edm::ParameterSet& iConfig):
  isMC_    (iConfig.getUntrackedParameter<bool>("sampleLabel",  true)),
  TTbarMC_ (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0))
{
  //now do what ever initialization is needed
  genToken_      = consumes<edm::View<reco::GenParticle> >  (iConfig.getParameter<edm::InputTag>("genLabel"));
  muonToken_     = consumes<edm::View<cat::Muon> >          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_ = consumes<edm::View<cat::Electron> >      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_      = consumes<edm::View<cat::Jet> >           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_      = consumes<edm::View<cat::MET> >           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_       = consumes<int>                            (iConfig.getParameter<edm::InputTag>("pvLabel"));
  puWeight_      = consumes<float>                          (iConfig.getParameter<edm::InputTag>("puWeight"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

  b_Jet_px = new std::vector<float>;
  b_Jet_py = new std::vector<float>;
  b_Jet_pz = new std::vector<float>;
  b_Jet_E  = new std::vector<float>;

  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;

  b_Jet_smearedRes     = new std::vector<float>;
  b_Jet_smearedResDown = new std::vector<float>;
  b_Jet_smearedResUp   = new std::vector<float>;
  b_Jet_shiftedEnUp    = new std::vector<float>;
  b_Jet_shiftedEnDown  = new std::vector<float>;

  b_Jet_CSV  = new std::vector<float>;

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "Event/I");
  tree->Branch("run",        &b_Run,         "Run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "Lumi_Number/I");

  tree->Branch("PUWeight", &b_PUWeight, "PUWeight/F");
  tree->Branch("GoodPV",   &b_nGoodPV,  "nGoodPV/I");

  tree->Branch("channel",  &b_Channel,  "Channel/I");

  tree->Branch("MET",     &b_MET,     "MET/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("lepton_px", &b_Lepton_px, "lepton_px/F");
  tree->Branch("lepton_py", &b_Lepton_py, "lepton_py/F");
  tree->Branch("lepton_pz", &b_Lepton_pz, "lepton_pz/F");
  tree->Branch("lepton_E" , &b_Lepton_E,  "lepton_E/F" );

  tree->Branch("jet_px", "std::vector<float>", &b_Jet_px);
  tree->Branch("jet_py", "std::vector<float>", &b_Jet_py);
  tree->Branch("jet_pz", "std::vector<float>", &b_Jet_pz);
  tree->Branch("jet_E" , "std::vector<float>", &b_Jet_E );

  tree->Branch("jet_partonFlavour", "std::vector<int>", &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>", &b_Jet_hadronFlavour);

  tree->Branch("jet_smearedRes",     "std::vector<float>", &b_Jet_smearedRes);
  tree->Branch("jet_smearedResDown", "std::vector<float>", &b_Jet_smearedResDown);
  tree->Branch("jet_smearedResUp",   "std::vector<float>", &b_Jet_smearedResUp);
  tree->Branch("jet_shiftedEnUp",    "std::vector<float>", &b_Jet_shiftedEnUp);
  tree->Branch("jet_shiftedEnDown",  "std::vector<float>", &b_Jet_shiftedEnDown);

  tree->Branch("jet_CSV" , "std::vector<float>", &b_Jet_CSV );

}


TtbarSingleLeptonAnalyzer::~TtbarSingleLeptonAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete b_Jet_px;
  delete b_Jet_py;
  delete b_Jet_pz;
  delete b_Jet_E;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_smearedRes;
  delete b_Jet_smearedResDown;
  delete b_Jet_smearedResUp;
  delete b_Jet_shiftedEnUp;
  delete b_Jet_shiftedEnDown;

  delete b_Jet_CSV;

}

//
// member functions
//

// ------------ method called for each event  ------------
void TtbarSingleLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_Jet_px->clear();
  b_Jet_py->clear();
  b_Jet_pz->clear();
  b_Jet_E->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_smearedRes->clear();
  b_Jet_smearedResDown->clear();
  b_Jet_smearedResUp->clear();
  b_Jet_shiftedEnUp->clear();
  b_Jet_shiftedEnDown->clear();

  b_Jet_CSV->clear();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // PU Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  if(isMC_) {
    edm::Handle<float> PUWeight;

    iEvent.getByToken(puWeight_, PUWeight);
    b_PUWeight = *PUWeight;
  }
  else b_PUWeight = 1;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Generated Particles (For Pythia8)
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  int nGenLep = 999;
  bool GenLep_m = false;
  bool GenLep_p = false;

  if(TTbarMC_ > 0) {

    edm::Handle<edm::View<reco::GenParticle> > genParticles;
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
	
          if (id == -11 || id == -13) GenLep_m = true;
          if (id ==  11 || id ==  13) GenLep_p = true;
	
          // Taus
          if(abs(id) == 15){
            for(unsigned int h = 0; h <  gp.numberOfDaughters(); h++) {
              const reco::Candidate *gd = gp.daughter(h);
              const int taudauid = gd->pdgId();
              if (taudauid == -11 || taudauid == -13) GenLep_m = true;
              if (taudauid == 11 || taudauid == 13) GenLep_p = true;
            } // for(taus' daughters)
          } // if(taus)

        } // if(t->W)
	
      }// if (mu || e || tau)
    } // for(genParticles)

    if(!GenLep_p && !GenLep_m) nGenLep = 0; // Full Hadronic
    if((GenLep_p && !GenLep_m) || (!GenLep_p && GenLep_m)) nGenLep = 1; // Single Lepton
    else if(GenLep_p && GenLep_m) nGenLep = 2; // Dilepton
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

  Handle<edm::View<cat::MET> > MET;
  iEvent.getByToken(metToken_, MET);

  // MET-PF
  b_MET     = MET->begin()->pt();
  b_MET_phi = MET->begin()->phi();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  std::vector<cat::Electron> selectedElectrons;
  std::vector<cat::Electron> vetoElectrons;

  Handle<edm::View<cat::Electron> > electrons;
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

  std::vector<cat::Muon> selectedMuons;
  std::vector<cat::Muon> vetoMuons;

  Handle<edm::View<cat::Muon> > muons;
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

  bool EvTrigger = true; // Trigger requirement not yet applied
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  if ( (ch_tag == 0 && trigHelper.triggerFired("HLT_Mu24_eta2p1_v1")) ||
       (ch_tag == 1 && trigHelper.triggerFired("HLT_Ele27_eta2p1_WP75_Gsf_v1")) ) {
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

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Jets
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    Handle<edm::View<cat::Jet> > jets;
    iEvent.getByToken(jetToken_, jets);

    for (unsigned int i = 0; i < jets->size() ; i++) {

      const cat::Jet & jet = jets->at(i);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection
      if(std::abs(jet.eta()) < 2.4 && jet.pt() > 20 && jet.LooseId()) goodJet = true;
      // Jet Cleaning
      TLorentzVector vjet(jet.px(), jet.py(), jet.pz(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if(dr_LepJet > 0.4) cleanJet = true; 	

      if(goodJet && cleanJet){
	// Basic variables
	b_Jet_px->push_back(jet.px());
	b_Jet_py->push_back(jet.py());
	b_Jet_pz->push_back(jet.pz());
	b_Jet_E ->push_back(jet.energy());
	
	// Parton Flavour
	b_Jet_partonFlavour->push_back(jet.partonFlavour());
	b_Jet_hadronFlavour->push_back(jet.hadronFlavour());
	
	// Smeared and Shifted
	b_Jet_smearedRes     ->push_back(jet.smearedRes() );
	b_Jet_smearedResDown ->push_back(jet.smearedResDown());
	b_Jet_smearedResUp   ->push_back(jet.smearedResUp());
	b_Jet_shiftedEnUp    ->push_back(jet.shiftedEnUp());
	b_Jet_shiftedEnDown  ->push_back(jet.shiftedEnDown());
	
	// b-tag discriminant
	float jet_btagDis_CSV = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	b_Jet_CSV ->push_back(jet_btagDis_CSV);

	// BTagSFUtil *fBTagSF;   //The BTag SF utility
	// fBTagSF = new BTagSFUtil("CSVv2", "Medium", 0);

	// if (fBTagSF->IsTagged(jet_btagDis_CSV, -999999, jet.pt(), jet.eta()) ){
	//    std::cout << "Is btag Jet....." << std::endl;
	// }
      }
    }

    // Fill Tree with event at 1 lepton cut level
    tree->Fill();

  } // if(ch_tag)

}

//------------- Good Muon Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsTightMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Tight selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.isTightMuon());

  GoodMuon &= (i_muon_candidate.isPFMuon());       // PF
  GoodMuon &= (i_muon_candidate.pt()> 20);         // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.1); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.12 );

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodMuon;
}
//------------- Loose Muon Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsLooseMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Loose selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.isLooseMuon());

  GoodMuon &= (i_muon_candidate.isPFMuon());       // PF
  GoodMuon &= (i_muon_candidate.pt()> 15);         // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.12 );

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsTightElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );            // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);          // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.1);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
		   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose") > 0.0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;

  GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.12 );

  // Effective Area Parametrization can be found in:
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}
//------------- Loose Electron Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsLooseElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );            // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);          // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
		   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0.0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.12 );
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarSingleLeptonAnalyzer);
