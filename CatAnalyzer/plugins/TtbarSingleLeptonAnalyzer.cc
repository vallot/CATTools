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
// Original Author:  Javier Brochero Cifuentes,512 1-001,+41227670488,
//         Created:  Tue Feb  3 09:52:55 CET 2015
// $Id$
//
//

// system include files
#include <memory>
#include <math.h> 
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm> // max

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
//
// class declaration
//

class TtbarSingleLeptonAnalyzer : public edm::EDAnalyzer {
public:
  explicit TtbarSingleLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarSingleLeptonAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  bool IsTightMuon(const cat::Muon & i_muon_candidate);
  bool IsTightElectron(const cat::Electron & i_electron_candidate);
  
  edm::EDGetTokenT<edm::View<cat::Muon> > muonToken_;
  edm::EDGetTokenT<edm::View<cat::Electron> > electronToken_;
  edm::EDGetTokenT<edm::View<cat::Jet> > jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> > metToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > pvToken_;
  edm::EDGetTokenT<double> puWeight_;

  // ----------member data ---------------------------

  TTree *vallot = new TTree();

  unsigned int minTracks_;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Tree Branches
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Event info
  int b_Event, b_Run, b_Lumi_Number;

  // PU/Vertices
  float b_PUWeight; 
  int b_nPV, b_nGoodPV;

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
  // ID
  std::vector<bool> *b_Jet_LooseID;
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
TtbarSingleLeptonAnalyzer::TtbarSingleLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonToken_ = consumes<edm::View<cat::Muon> >(iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_ = consumes<edm::View<cat::Electron> >(iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_ = consumes<edm::View<cat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_ = consumes<edm::View<cat::MET> >(iConfig.getParameter<edm::InputTag>("metLabel"));     
  pvToken_ = consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("pvLabel"));
  puWeight_ = consumes<double>(iConfig.getParameter<edm::InputTag>("puWeight"));

  edm::Service<TFileService> fs;
  vallot = fs->make<TTree>("vallot", "TopTree");
  
  vallot->Branch("event",      &b_Event,       "Event/I");
  vallot->Branch("run",        &b_Run,         "Run/I");
  vallot->Branch("luminumber", &b_Lumi_Number, "Lumi_Number/I");

  vallot->Branch("PUWeight", &b_PUWeight, "PUWeight/F");
  vallot->Branch("PV",       &b_nPV,      "nPV/I");
  vallot->Branch("GoodPV",   &b_nGoodPV,  "nGoodPV/I");

  vallot->Branch("channel",  &b_Channel,  "Channel/I");

  vallot->Branch("MET",     &b_MET,     "MET/F");
  vallot->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  vallot->Branch("lepton_px", &b_Lepton_px, "lepton_px/F");
  vallot->Branch("lepton_py", &b_Lepton_py, "lepton_py/F");
  vallot->Branch("lepton_pz", &b_Lepton_pz, "lepton_pz/F");
  vallot->Branch("lepton_E" , &b_Lepton_E,  "lepton_E/F" );

  vallot->Branch("jet_px", "std::vector<float>", &b_Jet_px);
  vallot->Branch("jet_py", "std::vector<float>", &b_Jet_py);
  vallot->Branch("jet_pz", "std::vector<float>", &b_Jet_pz);
  vallot->Branch("jet_E" , "std::vector<float>", &b_Jet_E );

  vallot->Branch("jet_LooseID", "std::vector<bool>", &b_Jet_LooseID);
  
  vallot->Branch("jet_partonFlavour", "std::vector<int>", &b_Jet_partonFlavour);
  vallot->Branch("jet_hadronFlavour", "std::vector<int>", &b_Jet_hadronFlavour);
  
  vallot->Branch("jet_smearedRes",     "std::vector<float>", &b_Jet_smearedRes);
  vallot->Branch("jet_smearedResDown", "std::vector<float>", &b_Jet_smearedResDown);
  vallot->Branch("jet_smearedResUp",   "std::vector<float>", &b_Jet_smearedResUp); 
  vallot->Branch("jet_shiftedEnUp",    "std::vector<float>", &b_Jet_shiftedEnUp);  
  vallot->Branch("jet_shiftedEnDown",  "std::vector<float>", &b_Jet_shiftedEnDown);

  vallot->Branch("jet_CSV" , "std::vector<float>", &b_Jet_CSV );

}


TtbarSingleLeptonAnalyzer::~TtbarSingleLeptonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void TtbarSingleLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  b_Jet_px = new std::vector<float>;
  b_Jet_py = new std::vector<float>;
  b_Jet_pz = new std::vector<float>;
  b_Jet_E  = new std::vector<float>;

  b_Jet_LooseID = new std::vector<bool>;
  
  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;
  
  b_Jet_smearedRes     = new std::vector<float>;
  b_Jet_smearedResDown = new std::vector<float>;
  b_Jet_smearedResUp   = new std::vector<float>;
  b_Jet_shiftedEnUp    = new std::vector<float>;
  b_Jet_shiftedEnDown  = new std::vector<float>;
  
  b_Jet_CSV  = new std::vector<float>;
  

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
  edm::Handle<double> PUWeight;

  iEvent.getByToken(puWeight_, PUWeight);
  b_PUWeight = *PUWeight;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Primary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  edm::Handle< std::vector<reco::Vertex> > pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );
  reco::Vertex primaryVertex;

  int n_vtxs = 0;
  // Loop over vertices
  if (pvHandle->size() != 0) {
    for (unsigned int i=0; i< pvHandle->size(); i++) {
     reco::Vertex vtx = pvHandle->at(i);
     if ( fabs(vtx.z())        < 24 &&
	  vtx.position().Rho() < 2  &&
	  vtx.ndof()           > 4  &&
	  !(vtx.isFake())              ) {
       n_vtxs ++;
     }
    }// for(vertex)
  } // if(vertex)
  
  b_nPV = pvHandle->size();
  b_nGoodPV = n_vtxs;  

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Secondary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  
  //do we need this?
 
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

  Handle<edm::View<cat::Electron> > electrons;
  iEvent.getByToken(electronToken_, electrons); 
 
  for (unsigned int i = 0; i < electrons->size() ; i++) { 
    const cat::Electron & electron = electrons->at(i);
    if( IsTightElectron( electron ) ) selectedElectrons.push_back( electron );
  }


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  std::vector<cat::Muon> selectedMuons; 

  Handle<edm::View<cat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons); 
 
  for (unsigned int i = 0; i < muons->size() ; i++) { 
    const cat::Muon & muon = muons->at(i);
    if( IsTightMuon( muon) ) selectedMuons.push_back( muon);
  }

  //---------------------------------------------------------------------------
  //----------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  TLorentzVector lepton;
  int ch_tag  =999;

  if(selectedMuons.size() == 1 && selectedElectrons.size() == 0){
    lepton.SetPxPyPzE(selectedMuons[0].px(), selectedMuons[0].py(), selectedMuons[0].pz(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets
  }

  if(selectedMuons.size() == 0 && selectedElectrons.size() == 1){
    lepton.SetPxPyPzE(selectedElectrons[0].px(), selectedElectrons[0].py(), selectedElectrons[0].pz(), selectedElectrons[0].energy());
    ch_tag = 1; //electron + jets
  }


  //---------------------------------------------------------------------------
  //----------------------------------------------------------------
  // Fill Tree with events that have ONLY one lepton
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  if (ch_tag<2){ // Single lepton event 

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

        if(fabs(jet.eta()) < 2.4 && jet.pt() > 20 ){
	
	    // Basic variables
	    b_Jet_px->push_back(jet.px());
	    b_Jet_py->push_back(jet.py());
	    b_Jet_pz->push_back(jet.pz());
	    b_Jet_E ->push_back(jet.energy());

	    // Jet ID (Loose)
	    b_Jet_LooseID ->push_back(jet.LooseId());

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
	    float jet_btagDis_CSV = jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
	    b_Jet_CSV ->push_back(jet_btagDis_CSV);
	
        }
    }
    
    vallot->Fill();
    
  } // if(ch_tag)

  delete b_Jet_px;
  delete b_Jet_py;
  delete b_Jet_pz;
  delete b_Jet_E;

  delete b_Jet_LooseID;
  
  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;
  
  delete b_Jet_smearedRes;
  delete b_Jet_smearedResDown;
  delete b_Jet_smearedResUp;
  delete b_Jet_shiftedEnUp;
  delete b_Jet_shiftedEnDown;

  delete b_Jet_CSV;
 

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByToken("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

//------------- Good Muon Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsTightMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;
  
  // Tight cut already defined into CAT::Muon
  //GoodMuon &= (i_muon_candidate->isTightMuon());
  
  GoodMuon &= (i_muon_candidate.isPFMuon());       // PF
  GoodMuon &= (i_muon_candidate.pt()> 20);         // pT
  GoodMuon &= (fabs(i_muon_candidate.eta())< 2.4); // eta

  GoodMuon &=(i_muon_candidate.isGlobalMuon());
  GoodMuon &=(i_muon_candidate.isPFMuon());  
  GoodMuon &=(i_muon_candidate.normalizedChi2() < 10);  
  GoodMuon &=(i_muon_candidate.numberOfValidMuonHits() > 0);  
  GoodMuon &=(i_muon_candidate.numberOfMatchedStations() > 1);  
  GoodMuon &=(fabs(i_muon_candidate.dxy()) < 0.2); //mm
  GoodMuon &=(fabs(i_muon_candidate.dz()) < 0.5); //mm
  GoodMuon &=(i_muon_candidate.numberOfValidPixelHits() > 0);
  GoodMuon &=(i_muon_candidate.trackerLayersWithMeasurement() > 5);

  float PFIsoMuon=999.;
  PFIsoMuon = i_muon_candidate.chargedHadronIso(0.3) +
              std::max(0.0, i_muon_candidate.neutralHadronIso(0.3) + 
	                    i_muon_candidate.photonIso(0.3) - 
	                    0.5*i_muon_candidate.puChargedHadronIso(0.3));
    
  PFIsoMuon = PFIsoMuon/i_muon_candidate.pt();
  
  GoodMuon &=( PFIsoMuon<0.12 );

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsTightElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );            // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);          // pT
  GoodElectron &= (fabs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (fabs(i_electron_candidate.eta()) < 1.4442 || 
		   fabs(i_electron_candidate.eta()) > 1.566);

  GoodElectron &= i_electron_candidate.passConversionVeto();

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium") > 0.0;

//----------------------------------------------------------------------------------------------------
//------------- The Relative Isolation is already calculated in the CAT object -----------------------
//----------------------------------------------------------------------------------------------------
  // Effective Area Parametrization
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  // Double_t AEff03 = 0.;

  // if      (fabs(i_electron_candidate->eta()) < 0.8)                                                 AEff03 = 0.1013;
  // else if (fabs(i_electron_candidate->eta()) >= 0.8   && fabs(i_electron_candidate->eta()) < 1.3)   AEff03 = 0.0988; 
  // else if (fabs(i_electron_candidate->eta()) >= 1.3   && fabs(i_electron_candidate->eta()) < 2.0)   AEff03 = 0.0572; 
  // else if (fabs(i_electron_candidate->eta()) >= 2.0   && fabs(i_electron_candidate->eta()) < 2.2)   AEff03 = 0.0842; 
  // else if (fabs(i_electron_candidate->eta()) >= 2.2)                                                AEff03 = 0.1530; 
  
  // float PFIsoElectron = ( i_electron_candidate->chargedHadronIso( 0.3 ) +
  // 			  std::max(0.0, 
  // 				   i_electron_candidate->neutralHadronIso( 0.3 ) +
  // 				   i_electron_candidate->photonIso( 0.3 ) -  
  // 				   AEff03*1
  // 				   )
  // 			  );
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------


  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;
  GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.12 );

  return GoodElectron;

}
// ------------ method called once each job just before starting event loop  ------------
void 
TtbarSingleLeptonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TtbarSingleLeptonAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TtbarSingleLeptonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TtbarSingleLeptonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TtbarSingleLeptonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TtbarSingleLeptonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TtbarSingleLeptonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarSingleLeptonAnalyzer);
