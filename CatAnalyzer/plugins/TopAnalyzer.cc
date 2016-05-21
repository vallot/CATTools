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
#include "CATTools/DataFormats/interface/GenWeights.h"

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
      edm::EDGetTokenT<cat::GenWeights>                genWeightToken_;

      // ----------member data ---------------------------

      TTree * tree;
      TH1F * tmp;
 
      int EVENT;
      int RUN;
      int LUMI;
    
      float PUWeight;
      float GenWeight;
      int NVertex; 

      double MET;
      double MET_Px;
      double MET_Py;

      const int kMax = 100;

      int NMuon;

      float Muon_Pt[100];
      float Muon_Eta[100];
      float Muon_Phi[100];
      float Muon_E[100];
      float Muon_Iso03[100];
      float Muon_Iso04[100];
      float Muon_Charge[100];

      int NLooseMuon;

      float LooseMuon_Pt[100];
      float LooseMuon_Eta[100];
      float LooseMuon_Phi[100];
      float LooseMuon_E[100];
      float LooseMuon_Iso03[100];
      float LooseMuon_Iso04[100];
      float LooseMuon_Charge[100];

      int NElectron;

      float Electron_Pt[100];
      float Electron_Eta[100];
      float Electron_Phi[100];
      float Electron_E[100];
      float Electron_Iso03[100];
      float Electron_Iso04[100];
      float Electron_Charge[100];

      int NLooseElectron;

      float LooseElectron_Pt[100];
      float LooseElectron_Eta[100];
      float LooseElectron_Phi[100];
      float LooseElectron_E[100];
      float LooseElectron_Iso03[100];
      float LooseElectron_Iso04[100];
      float LooseElectron_Charge[100];


      int NJet;
     
      float Jet_Pt[100];
      float Jet_Eta[100];
      float Jet_Phi[100];
      float Jet_E[100];
      float Jet_partonFlavour[100];
      float Jet_hadronFlavour[100];
      float Jet_BTag[100];
      float Jet_bDiscriminator[100];
 
      float Jet_JES_Up[100];
      float Jet_JES_Dw[100];
 
      int NBJet;

      int DiLeptonic;
      int SemiLeptonic;
 
      int GenNJet20;
      int GenNBJet20;
      int GenNCJet20;
      int GenNAddJet20;
      int GenNAddBJet20;
      int GenNAddCJet20;
      
      float GenLepton1_Pt;
      float GenLepton1_Eta;
      float GenLepton2_Pt;
      float GenLepton2_Eta;
 
      float WMuon_MT[100];
      float WMuon_Phi[100];
      float WElectron_MT[100];
      float WElectron_Phi[100]; 
      
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
   genTopToken_      = consumes<cat::GenTopCollection>       (iConfig.getParameter<edm::InputTag>("genTopLabel"));
   genToken_      = consumes<reco::GenParticleCollection>    (iConfig.getParameter<edm::InputTag>("genLabel"));
   muonToken_     = consumes<cat::MuonCollection>            (iConfig.getParameter<edm::InputTag>("muonLabel"));
   electronToken_ = consumes<cat::ElectronCollection>        (iConfig.getParameter<edm::InputTag>("electronLabel"));
   jetToken_      = consumes<cat::JetCollection>             (iConfig.getParameter<edm::InputTag>("jetLabel"));
   metToken_      = consumes<cat::METCollection>             (iConfig.getParameter<edm::InputTag>("metLabel"));
   pvToken_       = consumes<int>                            (iConfig.getParameter<edm::InputTag>("pvLabel"));
   puWeight_      = consumes<float>                          (iConfig.getParameter<edm::InputTag>("puWeight"));
   genWeightToken_  = consumes<cat::GenWeights>              (iConfig.getParameter<edm::InputTag>("genWeightLabel"));


   usesResource("TFileService");

   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("events", "Tree for Top quark study");
   tmp = fs->make<TH1F>("EventSummary","EventSummary",2,0,2);
 
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

   if( genTops.isValid() ) {
     cat::GenTop catGenTop = genTops->at(0);
     DiLeptonic = catGenTop.diLeptonic(0);
     SemiLeptonic = catGenTop.semiLeptonic(0);

     GenNJet20 = catGenTop.NJets20();
     GenNBJet20 = catGenTop.NbJets20();
     GenNCJet20 = catGenTop.NcJets20();
     GenNAddJet20 = catGenTop.NaddJets20();
     GenNAddBJet20 = catGenTop.NaddbJets20();
     GenNAddCJet20 = catGenTop.NaddcJets20();

     if( catGenTop.lepton1().pt() > 0){
       GenLepton1_Pt = catGenTop.lepton1().pt();
       GenLepton1_Eta = catGenTop.lepton1().eta();
     }

     if( catGenTop.lepton2().pt() > 0){
       GenLepton2_Pt = catGenTop.lepton2().pt();
       GenLepton2_Eta = catGenTop.lepton2().eta();
     }

   }

   if(!iEvent.isRealData()) {
     edm::Handle<float> PileUpWeight;
     iEvent.getByToken(puWeight_, PileUpWeight);
     PUWeight = *PileUpWeight;

     edm::Handle<cat::GenWeights> genWeight;
     iEvent.getByToken(genWeightToken_, genWeight);
     GenWeight = genWeight->genWeight();
     tmp->Fill(1,GenWeight);
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

     bool passLooseMuon = muon.pt() > 15 && fabs(muon.eta()) < 2.4 && muon.isLooseMuon();  
     //bool passID = muon.pt() > 30 && fabs(muon.eta()) < 2.1 && muon.isTightMuon(); 

     if( !passLooseMuon ) continue;
   
     LooseMuon_Pt[nmuons] = muon.pt();
     LooseMuon_Eta[nmuons] = muon.eta();
     LooseMuon_Phi[nmuons] = muon.phi();
     LooseMuon_E[nmuons] = muon.energy();
     LooseMuon_Iso03[nmuons] = muon.relIso(0.3);
     LooseMuon_Iso04[nmuons] = muon.relIso(0.4);
     LooseMuon_Charge[nmuons] = muon.charge();
     nloosemuons++;

     bool passTightMuon = muon.pt() > 30 && fabs(muon.eta()) < 2.1 && muon.isTightMuon();
     //bool passIso = muon.relIso(0.4) < 0.15;

     if( !passTightMuon ) continue;

     Muon_Pt[nmuons] = muon.pt(); 
     Muon_Eta[nmuons] = muon.eta(); 
     Muon_Phi[nmuons] = muon.phi(); 
     Muon_E[nmuons] = muon.energy();
     Muon_Iso03[nmuons] = muon.relIso(0.3);
     Muon_Iso04[nmuons] = muon.relIso(0.4);
     Muon_Charge[nmuons] = muon.charge();

     WMuon_MT[nmuons] = transverseMass( muon.p4(), METHandle->begin()->p4() );
     WMuon_Phi[nmuons] = fabs(deltaPhi( muon.phi(), METHandle->begin()->p4().phi()));
     nmuons++;

   }

   NMuon = nmuons;
   NLooseMuon = nloosemuons;
   int nelectrons = 0;
   int nlooseelectrons = 0;
   for (unsigned int i = 0; i < electrons->size() ; i++) {
     const cat::Electron & electron = electrons->at(i);

     bool passLooseElectron = electron.pt() > 15 && fabs(electron.eta()) < 2.4 && electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto") > 0;
     //bool passID = electron.pt() > 30 && fabs(electron.eta()) < 2.1 && electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") > 0;

     if ( !passLooseElectron ) continue;
  
     LooseElectron_Pt[nelectrons] = electron.pt();
     LooseElectron_Eta[nelectrons] = electron.eta();
     LooseElectron_Phi[nelectrons] = electron.phi();
     LooseElectron_E[nelectrons] = electron.energy();
     LooseElectron_Iso03[nelectrons] = electron.relIso(0.3);
     LooseElectron_Iso04[nelectrons] = electron.relIso(0.4);
     LooseElectron_Charge[nelectrons] = electron.charge();
     nlooseelectrons++;

     bool passMediumElectron = electron.pt() > 30 && fabs(electron.eta()) < 2.1 && electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0;
     //bool passIso = electron.relIso() < 0.12;

     if ( !passMediumElectron ) continue;

     Electron_Pt[nelectrons] = electron.pt();
     Electron_Eta[nelectrons] = electron.eta();
     Electron_Phi[nelectrons] = electron.phi();
     Electron_E[nelectrons] = electron.energy();
     Electron_Iso03[nelectrons] = electron.relIso(0.3);
     Electron_Iso04[nelectrons] = electron.relIso(0.4);
     Electron_Charge[nelectrons] = electron.charge();

     WElectron_MT[nelectrons] = transverseMass( electron.p4(), METHandle->begin()->p4() );
     WElectron_Phi[nelectrons] = fabs(deltaPhi( electron.phi(), METHandle->begin()->p4().phi()));
     nelectrons++;

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

     for(int j = 0 ; j < NMuon ; j++){ 
       if( Muon_Iso03[j] < 0.15 ){
         TLorentzVector vlep(Muon_Pt[j], Muon_Eta[j], Muon_Phi[j], Muon_E[j]);
         dr = vjet.DeltaR(vlep);
       }
     }

     for(int j = 0 ; j < NElectron ; j++){
       if( Electron_Iso03[j] < 0.15 ){
         TLorentzVector vlep(Electron_Pt[j], Electron_Eta[j], Electron_Phi[j], Electron_E[j]);
         dr = vjet.DeltaR(vlep);
       }
     }

     if( dr < 0.4) continue;

     Jet_Pt[nJets] = jet.pt(); 
     Jet_Eta[nJets] = jet.eta();
     Jet_Phi[nJets] = jet.phi();
     Jet_E[nJets] = jet.energy();

     Jet_partonFlavour[nJets] = jet.partonFlavour();
     Jet_hadronFlavour[nJets] = jet.hadronFlavour();

     Jet_JES_Up[nJets] = jet.shiftedEnUp();
     Jet_JES_Dw[nJets] = jet.shiftedEnDown();

     double bDiscriminator = jet.bDiscriminator(BTAG_CSVv2);
     Jet_bDiscriminator[nJets] = bDiscriminator;
     if( bDiscriminator > WP_BTAG_CSVv2M) {
       nbJets++;
       Jet_BTag[nJets] = 1;
     }else{
       Jet_BTag[nJets] = 0;
     }

     nJets++;
   }

   NJet = nJets;
   NBJet = nbJets;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   if (NMuon > 0 || NElectron > 0 ) {
     tree->Fill();
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TopAnalyzer::beginJob()
{
   tree->Branch("EVENT",&EVENT,"EVENT/i");
   tree->Branch("RUN",&RUN,"RUN/i");
   tree->Branch("LUMI",&LUMI,"LUMI/i");
   tree->Branch("PUWeight",&PUWeight,"PUWeight/F");
   tree->Branch("GenWeight",&GenWeight,"GenWeight/F");
   tree->Branch("NVertex",&NVertex,"NVertex/i");

   tree->Branch("MET",&MET,"MET/d");
   tree->Branch("MET_Px",&MET_Px,"MET_Px/d");
   tree->Branch("MET_Py",&MET_Py,"MET_Py/d");

   tree->Branch("NMuon",&NMuon,"NMuon/I");
   tree->Branch("Muon_Pt",Muon_Pt,"Muon_Pt[NMuon]/F");
   tree->Branch("Muon_Eta",Muon_Eta,"Muon_Eta[NMuon]/F");
   tree->Branch("Muon_Phi",Muon_Phi,"Muon_Phi[NMuon]/F");
   tree->Branch("Muon_E",Muon_E,"Muon_E[NMuon]/F");
   tree->Branch("Muon_Iso03",Muon_Iso03,"Muon_Iso03[NMuon]/F");
   tree->Branch("Muon_Iso04",Muon_Iso04,"Muon_Iso04[NMuon]/F");
   tree->Branch("Muon_Charge",Muon_Charge,"Muon_Charge[NMuon]/F");

   tree->Branch("NLooseMuon",&NLooseMuon,"NLooseMuon/I");
   tree->Branch("LooseMuon_Pt",LooseMuon_Pt,"LooseMuon_Pt[NLooseMuon]/F");
   tree->Branch("LooseMuon_Eta",LooseMuon_Eta,"LooseMuon_Eta[NLooseMuon]/F");
   tree->Branch("LooseMuon_Phi",LooseMuon_Phi,"LooseMuon_Phi[NLooseMuon]/F");
   tree->Branch("LooseMuon_E",LooseMuon_E,"LooseMuon_E[NLooseMuon]/F");
   tree->Branch("LooseMuon_Iso03",LooseMuon_Iso03,"LooseMuon_Iso03[NLooseMuon]/F");
   tree->Branch("LooseMuon_Iso04",LooseMuon_Iso04,"LooseMuon_Iso04[NLooseMuon]/F");
   tree->Branch("LooseMuon_Charge",LooseMuon_Charge,"LooseMuon_Charge[NLooseMuon]/F");

   tree->Branch("NElectron",&NElectron,"NElectron/I");
   tree->Branch("Electron_Pt",Electron_Pt,"Electron_Pt[NElectron]/F");
   tree->Branch("Electron_Eta",Electron_Eta,"Electron_Eta[NElectron]/F");
   tree->Branch("Electron_Phi",Electron_Phi,"Electron_Phi[NElectron]/F");
   tree->Branch("Electron_E",Electron_E,"Electron_E[NElectron]/F");
   tree->Branch("Electron_Iso03",Electron_Iso03,"Electron_Iso03[NElectron]/F");
   tree->Branch("Electron_Iso04",Electron_Iso04,"Electron_Iso04[NElectron]/F");
   tree->Branch("Electron_Charge",Electron_Charge,"Electron_Charge[NElectron]/F");

   tree->Branch("NLooseElectron",&NLooseElectron,"NLooseElectron/I");
   tree->Branch("LooseElectron_Pt",LooseElectron_Pt,"LooseElectron_Pt[NLooseElectron]/F");
   tree->Branch("LooseElectron_Eta",LooseElectron_Eta,"LooseElectron_Eta[NLooseElectron]/F");
   tree->Branch("LooseElectron_Phi",LooseElectron_Phi,"LooseElectron_Phi[NLooseElectron]/F");
   tree->Branch("LooseElectron_E",LooseElectron_E,"LooseElectron_E[NLooseElectron]/F");
   tree->Branch("LooseElectron_Iso03",LooseElectron_Iso03,"LooseElectron_Iso03[NLooseElectron]/F");
   tree->Branch("LooseElectron_Iso04",LooseElectron_Iso04,"LooseElectron_Iso04[NLooseElectron]/F");
   tree->Branch("LooseElectron_Charge",LooseElectron_Charge,"LooseElectron_Charge[NLooseElectron]/F");

   tree->Branch("NJet",&NJet,"NJet/i");
   tree->Branch("Jet_Pt",Jet_Pt,"Jet_Pt[NJet]/F");
   tree->Branch("Jet_Eta",Jet_Eta,"Jet_Eta[NJet]/F");
   tree->Branch("Jet_Phi",Jet_Phi,"Jet_Phi[NJet]/F");
   tree->Branch("Jet_E",Jet_E,"Jet_E[NJet]/F");
   tree->Branch("Jet_partonFlavour",Jet_partonFlavour,"Jet_partonFlavour[NJet]/F");
   tree->Branch("Jet_hadronFlavour",Jet_hadronFlavour,"Jet_hadronFlavour[NJet]/F");
   tree->Branch("Jet_BTag",Jet_BTag,"Jet_BTag[NJet]/F");
   tree->Branch("Jet_bDiscriminator",Jet_bDiscriminator,"Jet_bDiscriminator[NJet]/F"); 

   tree->Branch("Jet_JES_Up",Jet_JES_Up,"Jet_JES_Up[NJet]/F");
   tree->Branch("Jet_JES_Dw",Jet_JES_Dw,"Jet_JES_Dw[NJet]/F");
 
   tree->Branch("NBJet",&NBJet,"NBJet/i");
   tree->Branch("DiLeptonic",&DiLeptonic,"DiLeptonic/i");
   tree->Branch("SemiLeptonic",&SemiLeptonic,"SemiLeptonic/i");

   tree->Branch("GenNJet20",&GenNJet20, "GenNJet20/i");
   tree->Branch("GenNBJet20",&GenNBJet20, "GenNBJet20/i");
   tree->Branch("GenNCJet20",&GenNCJet20, "GenNCJet20/i");
   tree->Branch("GenNAddJet20",&GenNAddJet20, "GenNAddJet20/i");
   tree->Branch("GenNAddBJet20",&GenNAddBJet20, "GenNAddBJet20/i");
   tree->Branch("GenNAddCJet20",&GenNAddCJet20, "GenNAddCJet20/i");   

   tree->Branch("GenLepton1_Pt",&GenLepton1_Pt, "GenLepton1_Pt/f");
   tree->Branch("GenLepton1_Eta",&GenLepton1_Eta, "GenLepton1_Eta/f");
   tree->Branch("GenLepton2_Pt",&GenLepton2_Pt, "GenLepton2_Pt/f");
   tree->Branch("GenLepton2_Eta",&GenLepton2_Eta, "GenLepton2_Eta/f");

   tree->Branch("WMuon_MT",WMuon_MT,"WMuon_MT[NMuon]/F"); 
   tree->Branch("WMuon_Phi",WMuon_Phi,"WMuon_Phi[NMuon]/F"); 
   tree->Branch("WElectron_MT",WElectron_MT,"WElectron_MT[NElectron]/F"); 
   tree->Branch("WElectron_Phi",WElectron_Phi,"WElectron_Phi[NElectron]/F"); 
 
}

void
TopAnalyzer::clear(){

  PUWeight = 1.0;
  GenWeight = 1.0;
  NVertex = -1;
  NMuon = -1;
  NLooseMuon = -1;
  NElectron = -1;
  NLooseElectron = -1;
  NJet = -1;
  NBJet = -1;
  DiLeptonic = -1;
  SemiLeptonic = -1;
  
  GenNJet20 = -1; 
  GenNBJet20 = -1;
  GenNCJet20 = -1;
  GenNAddJet20 = -1;
  GenNAddBJet20 = -1;
  GenNAddCJet20 = -1;

  GenLepton1_Pt = -9.0;
  GenLepton1_Eta = -9.0;
  GenLepton2_Pt = -9.0;
  GenLepton2_Eta = -9.0;

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
