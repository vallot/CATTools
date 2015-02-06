// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

void analyze()
{
}

using namespace std;
//
// class decleration
//
class CATMuonAnalysis : public edm::EDAnalyzer {
  public:
    explicit CATMuonAnalysis(const edm::ParameterSet&);
    ~CATMuonAnalysis() {};

  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<edm::View<cat::Muon> > muonToken_;
    edm::EDGetTokenT<edm::View<cat::Electron> > electronToken_;
    edm::EDGetTokenT<edm::View<cat::Jet> > jetToken_;
    edm::EDGetTokenT<edm::View<cat::MET> > metToken_;

    TH1F* hMuon_phi_;
    TH1F* eta_;
    TH1F* pt_;
    TH1F* mass_;
    TH1F* chIso_;
    TH1F* nhIso_;
    TH1F* phIso_;
    TH1F* puIso_;
    
    TH1F* hElectron_phi;
};

CATMuonAnalysis::CATMuonAnalysis(const edm::ParameterSet& iConfig)
{
  muonToken_ = consumes<edm::View<cat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  electronToken_ = consumes<edm::View<cat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_ = consumes<edm::View<cat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_ = consumes<edm::View<cat::MET> >(iConfig.getParameter<edm::InputTag>("met"));

  edm::Service<TFileService> fs;
  TFileDirectory dirMuon = fs->mkdir("muon", "muon");
  hMuon_phi_ = dirMuon.make<TH1F>("phi","phi",100,-4,4);
  eta   = dirMuon.make<TH1F>("eta","eta",400,-3,3);
  pt    = dirMuon.make<TH1F>("pt","pt",500,0,500);
  mass  = dirMuon.make<TH1F>("mass","mass",1000,0,100);
  chIso = dirMuon.make<TH1F>("chIso","chIso",400,0,4);
  nhIso = dirMuon.make<TH1F>("nhIso","nhIso",400,0,4);
  phIso = dirMuon.make<TH1F>("phIso","phIso",400,0,4);
  puIso = dirMuon.make<TH1F>("puIso","puIso",400,0,4);

  TFileDirectory dirElectron = fs->mkdir("electron", "electron");
  hElectron_phi   = dirElectron.make<TH1F>("phi","phi",100,-4,4);
  eta   = dirElectron.make<TH1F>("eta","eta",400,-3,3);
  pt    = dirElectron.make<TH1F>("pt","pt",500,0,500);
  mass  = dirElectron.make<TH1F>("mass","mass",1000,0,100);
  chIso = dirElectron.make<TH1F>("chIso","chIso",400,0,4);
  nhIso = dirElectron.make<TH1F>("nhIso","nhIso",400,0,4);
  phIso = dirElectron.make<TH1F>("phIso","phIso",400,0,4);
  puIso = dirElectron.make<TH1F>("puIso","puIso",400,0,4);
}

void CATMuonAnalysis::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  edm::Handle<edm::View<cat::Muon> > muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  // Do the remainign things for electrons, jets, met
  // for dong : learn how to use pointers and iterators
  //for ( const cat::Muon& muon : *muonHandle )  // For c++0x
  //for ( edm::View<cat::Muon>::const_iterator iMuon = muonHandle->begin();
  //      iMuon != muonHandle->end(); ++iMuon )
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    //const cat::Muon& muon = *iMuon;
    const cat::Muon& muon = muonHandle->at(i);

    hMuon_phi_->Fill( muon.phi() );
    eta  ->Fill( muon.eta() );
    pt   ->Fill( muon.pt() );
    mass ->Fill( muon.mass() ); 
    chIso->Fill( muon.chargedHadronIso() );
    nhIso->Fill( muon.neutralHadronIso() );
    phIso->Fill( muon.photonIso() );
    puIso->Fill( muon.puChargedHadronIso() );
  }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATMuonAnalysis);

