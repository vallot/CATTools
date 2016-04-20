// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "TH1.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TopAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TopAnalyzer(const edm::ParameterSet&);
      ~TopAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      void clear();

      double transverseMass( const reco::Candidate::LorentzVector& lepton, const reco::Candidate::LorentzVector& met);

      edm::EDGetTokenT<cat::GenTopCollection>          genTopToken_;
      edm::EDGetTokenT<reco::GenParticleCollection>    genToken_;
      edm::EDGetTokenT<cat::MuonCollection>            muonToken_;
      edm::EDGetTokenT<cat::ElectronCollection>        electronToken_;
      edm::EDGetTokenT<cat::JetCollection>             jetToken_;
      edm::EDGetTokenT<cat::METCollection>             metToken_;
      edm::EDGetTokenT<int>                            pvToken_;
      edm::EDGetTokenT<float>                          puWeight_;

      // ----------member data ---------------------------

      TTree * tree;
      TH1F * tmp;

      int EVENT;
      int RUN;
      int LUMI;
    
      float PUWeight;
      int NVertex; 

      float MET;
      float MET_Px;
      float MET_Py;

      const int kMax = 100;

      int NMuon;

      std::vector<float> * Muon_Pt;
      std::vector<float> * Muon_Eta;
      std::vector<float> * Muon_Phi;
      std::vector<float> * Muon_E;
      std::vector<float> * Muon_Iso;
      std::vector<float> * Muon_Charge;

      int NLooseMuon;

      std::vector<float> * LooseMuon_Pt;
      std::vector<float> * LooseMuon_Eta;
      std::vector<float> * LooseMuon_Phi;
      std::vector<float> * LooseMuon_E;
      std::vector<float> * LooseMuon_Iso;
      std::vector<float> * LooseMuon_Charge;

      int NElectron;

      std::vector<float> * Electron_Pt;
      std::vector<float> * Electron_Eta;
      std::vector<float> * Electron_Phi;
      std::vector<float> * Electron_E;
      std::vector<float> * Electron_Iso;
      std::vector<float> * Electron_Charge;

      int NLooseElectron;

      std::vector<float> * LooseElectron_Pt;
      std::vector<float> * LooseElectron_Eta;
      std::vector<float> * LooseElectron_Phi;
      std::vector<float> * LooseElectron_E;
      std::vector<float> * LooseElectron_Iso;
      std::vector<float> * LooseElectron_Charge;


      int NJet;
     
      std::vector<float> * Jet_Pt;
      std::vector<float> * Jet_Eta;
      std::vector<float> * Jet_Phi;
      std::vector<float> * Jet_E;
      std::vector<float> * Jet_BTag;
      std::vector<float> * Jet_bDiscriminator;
  
      int NBJet;

      int DiLeptonic;
      int SemiLeptonic;
  
      std::vector<float> * WMuon_MT;
      std::vector<float> * WMuon_Phi;
      std::vector<float> * WElectron_MT;
      std::vector<float> * WElectron_Phi; 
      
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
TopAnalyzer::TopAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   genTopToken_      = consumes<cat::GenTopCollection>  (iConfig.getParameter<edm::InputTag>("genTopLabel"));
   genToken_      = consumes<reco::GenParticleCollection>  (iConfig.getParameter<edm::InputTag>("genLabel"));
   muonToken_     = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
   electronToken_ = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
   jetToken_      = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
   metToken_      = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
   pvToken_       = consumes<int>                            (iConfig.getParameter<edm::InputTag>("pvLabel"));
   puWeight_      = consumes<float>                          (iConfig.getParameter<edm::InputTag>("puWeight"));

   usesResource("TFileService");

   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("events", "Tree for Top quark study");
   tmp = fs->make<TH1F>("EventSummary","EventSummary",1,0,1);

   Muon_Pt = new std::vector<float>;
   Muon_Eta = new std::vector<float>;
   Muon_Phi = new std::vector<float>;
   Muon_E = new std::vector<float>;
   Muon_Iso = new std::vector<float>;
   Muon_Charge = new std::vector<float>;

   LooseMuon_Pt = new std::vector<float>;
   LooseMuon_Eta = new std::vector<float>;
   LooseMuon_Phi = new std::vector<float>;
   LooseMuon_E = new std::vector<float>;
   LooseMuon_Iso = new std::vector<float>;
   LooseMuon_Charge = new std::vector<float>;

   Electron_Pt = new std::vector<float>;
   Electron_Eta = new std::vector<float>;
   Electron_Phi = new std::vector<float>;
   Electron_E = new std::vector<float>;
   Electron_Iso = new std::vector<float>;
   Electron_Charge = new std::vector<float>;

   LooseElectron_Pt = new std::vector<float>;
   LooseElectron_Eta = new std::vector<float>;
   LooseElectron_Phi = new std::vector<float>;
   LooseElectron_E = new std::vector<float>;
   LooseElectron_Iso = new std::vector<float>;
   LooseElectron_Charge = new std::vector<float>;

   Jet_Pt = new std::vector<float>;
   Jet_Eta = new std::vector<float>;
   Jet_Phi = new std::vector<float>;
   Jet_E = new std::vector<float>; 
   Jet_BTag = new std::vector<float>; 
   Jet_bDiscriminator = new std::vector<float>; 

   WMuon_MT = new std::vector<float>;
   WMuon_Phi = new std::vector<float>;
   WElectron_MT = new std::vector<float>;
   WElectron_Phi = new std::vector<float>;
 
}


TopAnalyzer::~TopAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TopAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace cat;

   tmp->Fill(0); 

   clear();
 
   EVENT  = iEvent.id().event();
   RUN    = iEvent.id().run();
   LUMI   = iEvent.id().luminosityBlock();

   edm::Handle<cat::GenTopCollection> genTops;
   iEvent.getByToken(genTopToken_, genTops);

   cat::GenTop catGenTop = genTops->at(0);
   DiLeptonic = catGenTop.diLeptonic(0);
   SemiLeptonic = catGenTop.semiLeptonic(0);

   if(!iEvent.isRealData()) {
     edm::Handle<float> PileUpWeight;
     iEvent.getByToken(puWeight_, PileUpWeight);
     PUWeight = *PileUpWeight;
   }

   edm::Handle<int> pvHandle;
   iEvent.getByToken( pvToken_, pvHandle );

   NVertex = *pvHandle;

   Handle<cat::METCollection> METHandle;
   iEvent.getByToken(metToken_, METHandle);

   Handle<cat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
 
   Handle<cat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   Handle<cat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);

   MET     = METHandle->begin()->pt();
   MET_Px  = METHandle->begin()->px();
   MET_Py  = METHandle->begin()->py();

   int nmuons = 0;
   int nloosemuons = 0;
   for (unsigned int i = 0; i < muons->size() ; i++) {
     const cat::Muon & muon = muons->at(i);
     bool pass = muon.pt() > 30 && fabs(muon.eta()) < 2.1;
     if( !pass ) continue; 
     if( muon.isTightMuon() ){
       Muon_Pt->push_back(muon.pt()); 
       Muon_Eta->push_back(muon.eta()); 
       Muon_Phi->push_back(muon.phi()); 
       Muon_E->push_back(muon.energy());
       Muon_Iso->push_back(muon.relIso());
       Muon_Charge->push_back(muon.charge());

       WMuon_MT->push_back( transverseMass( muon.p4(), METHandle->begin()->p4() ) );
       WMuon_Phi->push_back( fabs( deltaPhi( muon.phi(), METHandle->begin()->p4().phi() )));
       nmuons++;
     }else if (muon.isLooseMuon()) {
       LooseMuon_Pt->push_back(muon.pt());
       LooseMuon_Eta->push_back(muon.eta());
       LooseMuon_Phi->push_back(muon.phi());
       LooseMuon_E->push_back(muon.energy());
       LooseMuon_Iso->push_back(muon.relIso());
       LooseMuon_Charge->push_back(muon.charge());
       nloosemuons++;
     }else{
       continue;
     }
   }
   NMuon = nmuons;
   NLooseMuon = nloosemuons;

   int nelectrons = 0;
   int nlooseelectrons = 0;
   for (unsigned int i = 0; i < electrons->size() ; i++) {
     const cat::Electron & electron = electrons->at(i);
     bool pass = electron.pt() > 30 && fabs(electron.eta()) < 2.1;
     if( !pass ) continue;
     if( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") > 0 ){
       Electron_Pt->push_back(electron.pt());
       Electron_Eta->push_back(electron.eta());
       Electron_Phi->push_back(electron.phi());
       Electron_E->push_back(electron.energy());
       Electron_Iso->push_back(electron.relIso());
       Electron_Charge->push_back(electron.charge());

       WElectron_MT->push_back(transverseMass( electron.p4(), METHandle->begin()->p4() ));
       WElectron_Phi->push_back(fabs(deltaPhi( electron.phi(), METHandle->begin()->p4().phi())));
       nelectrons++;
     }else if( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto") > 0  ) {
       LooseElectron_Pt->push_back(electron.pt());
       LooseElectron_Eta->push_back(electron.eta());
       LooseElectron_Phi->push_back(electron.phi());
       LooseElectron_E->push_back(electron.energy());
       LooseElectron_Iso->push_back(electron.relIso());
       LooseElectron_Charge->push_back(electron.charge());
       nlooseelectrons++;
     }
   }

   NElectron = nelectrons;
   NLooseElectron = nlooseelectrons;

   int nJets = 0;
   int nbJets = 0;

   for (unsigned int i = 0; i < jets->size() ; i++) {

     const cat::Jet & jet = jets->at(i);

     bool pass = std::abs(jet.eta()) < 2.4 && jet.pt() > 30 && jet.LooseId() ;
     if (!pass ) continue; 

     double dr = 999.9;
     TLorentzVector vjet(jet.px(), jet.py(), jet.pz(), jet.energy());

     double iso_cut = 0.12;

     for(int j = 0 ; j < NMuon ; j++){ 
       if( Muon_Iso->at(j) < iso_cut ){
         TLorentzVector vlep(Muon_Pt->at(j), Muon_Eta->at(j), Muon_Phi->at(j), Muon_E->at(j));
         dr = vjet.DeltaR(vlep);
       }
     }

     for(int j = 0 ; j < NElectron ; j++){
       if( Electron_Iso->at(j) < iso_cut ){
         TLorentzVector vlep(Electron_Pt->at(j), Electron_Eta->at(j), Electron_Phi->at(j), Electron_E->at(j));
         dr = vjet.DeltaR(vlep);
       }
     }

     if( dr < 0.4) continue;
     nJets++;

     Jet_Pt->push_back(jet.pt()); 
     Jet_Eta->push_back(jet.eta());
     Jet_Phi->push_back(jet.phi());
     Jet_E->push_back(jet.energy());

     double bDiscriminator = jet.bDiscriminator(BTAG_CSVv2);
     Jet_bDiscriminator->push_back(bDiscriminator);
     if( bDiscriminator > WP_BTAG_CSVv2M) {
       nbJets++;
       Jet_BTag->push_back(1);
     }else{
       Jet_BTag->push_back(0);
     }
   }

   NJet = nJets;
   NBJet = nbJets;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif


   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
TopAnalyzer::beginJob()
{
   tree->Branch("EVENT",&EVENT,"EVENT/i");
   tree->Branch("RUN",&RUN,"RUN/i");
   tree->Branch("LUMI",&LUMI,"LUMI/i");
   tree->Branch("PUWeight",&PUWeight,"PUWeight/i");
   tree->Branch("NVertex",&NVertex,"NVertex/i");

   tree->Branch("MET",&MET,"MET/d");
   tree->Branch("MET_Px",&MET_Px,"MET_Px/d");
   tree->Branch("MET_Py",&MET_Py,"MET_Py/d");

   tree->Branch("NMuon",&NMuon,"NMuon/I");
   tree->Branch("Muon_Pt","std::vector<float>",&Muon_Pt);
   tree->Branch("Muon_Eta","std::vector<float>",&Muon_Eta);
   tree->Branch("Muon_Phi","std::vector<float>",&Muon_Phi);
   tree->Branch("Muon_E","std::vector<float>",&Muon_E);
   tree->Branch("Muon_Iso","std::vector<float>",&Muon_Iso);
   tree->Branch("Muon_Chage","std::vector<float>",&Muon_Charge);

   tree->Branch("NLooseMuon",&NLooseMuon,"NLooseMuon/I");
   tree->Branch("LooseMuon_Pt","std::vector<float>",&LooseMuon_Pt);
   tree->Branch("LooseMuon_Eta","std::vector<float>",&LooseMuon_Eta);
   tree->Branch("LooseMuon_Phi","std::vector<float>",&LooseMuon_Phi);
   tree->Branch("LooseMuon_E","std::vector<float>",&LooseMuon_E);
   tree->Branch("LooseMuon_Iso","std::vector<float>",&LooseMuon_Iso);
   tree->Branch("LooseMuon_Chage","std::vector<float>",&LooseMuon_Charge);

   tree->Branch("NElectron",&NElectron,"NElectron/I");
   tree->Branch("Electron_Pt","std::vector<float>",&Electron_Pt);
   tree->Branch("Electron_Eta","std::vector<float>",&Electron_Eta);
   tree->Branch("Electron_Phi","std::vector<float>",&Electron_Phi);
   tree->Branch("Electron_E","std::vector<float>",&Electron_E);
   tree->Branch("Electron_Iso","std::vector<float>",&Electron_Iso);
   tree->Branch("Electron_Chage","std::vector<float>",&Electron_Charge);

   tree->Branch("NLooseElectron",&NLooseElectron,"NLooseElectron/I");
   tree->Branch("LooseElectron_Pt","std::vector<float>",&LooseElectron_Pt);
   tree->Branch("LooseElectron_Eta","std::vector<float>",&LooseElectron_Eta);
   tree->Branch("LooseElectron_Phi","std::vector<float>",&LooseElectron_Phi);
   tree->Branch("LooseElectron_E","std::vector<float>",&LooseElectron_E);
   tree->Branch("LooseElectron_Iso","std::vector<float>",&LooseElectron_Iso);
   tree->Branch("LooseElectron_Chage","std::vector<float>",&LooseElectron_Charge);

   tree->Branch("NJet",&NJet,"NJet/i");
   tree->Branch("Jet_Pt","std::vector<float>",&Jet_Pt);
   tree->Branch("Jet_Eta","std::vector<float>",&Jet_Eta);
   tree->Branch("Jet_Phi","std::vector<float>",&Jet_Phi);
   tree->Branch("Jet_E","std::vector<float>",&Jet_E);
   tree->Branch("Jet_BTag","std::vector<float>",&Jet_BTag); 
   tree->Branch("Jet_bDiscriminator","std::vector<float>",&Jet_bDiscriminator);  
 
   tree->Branch("NBJet",&NBJet,"NBJet/i");
   tree->Branch("DiLeptonic",&DiLeptonic,"DiLeptonic/i");
   tree->Branch("SemiLeptonic",&SemiLeptonic,"SemiLeptonic/i");

   tree->Branch("WMuon_MT","std::vector<float>",&WMuon_MT); 
   tree->Branch("WMuon_Phi","std::vector<float>",&WMuon_Phi); 
   tree->Branch("WElectron_MT","std::vector<float>",&WElectron_MT); 
   tree->Branch("WElectron_Phi","std::vector<float>",&WElectron_Phi); 
 

}

void
TopAnalyzer::clear(){

  PUWeight = 1.0;
  NVertex = -1;
  NMuon = -1;
  NLooseMuon = -1;
  NElectron = -1;
  NLooseElectron = -1;
  NJet = -1;
  NBJet = -1;
  DiLeptonic = -1;
  SemiLeptonic = -1;

  Muon_Pt->clear();
  Muon_Eta->clear();
  Muon_Phi->clear();
  Muon_E->clear();
  Muon_Iso->clear();
  Muon_Charge->clear();

  LooseMuon_Pt->clear();
  LooseMuon_Eta->clear();
  LooseMuon_Phi->clear();
  LooseMuon_E->clear();
  LooseMuon_Iso->clear();
  LooseMuon_Charge->clear();

  Electron_Pt->clear();
  Electron_Eta->clear();
  Electron_Phi->clear();
  Electron_E->clear();
  Electron_Iso->clear();
  Electron_Charge->clear();

  LooseElectron_Pt->clear();
  LooseElectron_Eta->clear();
  LooseElectron_Phi->clear();
  LooseElectron_E->clear();
  LooseElectron_Iso->clear();
  LooseElectron_Charge->clear();

  Jet_Pt->clear();
  Jet_Eta->clear();
  Jet_Phi->clear();
  Jet_E->clear(); 
  Jet_BTag->clear();
  Jet_bDiscriminator->clear();

  WMuon_MT->clear();
  WMuon_Phi->clear();
  WElectron_MT->clear();
  WElectron_Phi->clear();

}

double TopAnalyzer::transverseMass( const reco::Candidate::LorentzVector& lepton,
                                     const reco::Candidate::LorentzVector& met) {
  reco::Candidate::LorentzVector leptonT(lepton.Px(),lepton.Py(),0.,lepton.E()*sin(lepton.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+met;
  return std::sqrt(sumT.M2());
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TopAnalyzer::endJob() 
{


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TopAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TopAnalyzer);
