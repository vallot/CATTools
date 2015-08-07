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
#include "DataFormats/Common/interface/View.h"

#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Photon.h"
//#include "CATTools/DataFormats/interface/SecVertex.h"
#include "CATTools/DataFormats/interface/Tau.h"

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
class CATHisAnalysis : public edm::EDAnalyzer {
  public:
    explicit CATHisAnalysis(const edm::ParameterSet & pset);
    ~CATHisAnalysis() {};
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:

    // ----------member data ---------------------------
    typedef cat::Electron TElectron;
    typedef cat::GenJet TGenJet;
    typedef cat::Jet TJet;
    typedef cat::MET TMET;
    typedef cat::Muon TMuon;
    typedef cat::Photon TPhoton;
    typedef cat::Tau TTau;

    edm::EDGetTokenT<edm::View<TElectron> > electronToken_;
    edm::EDGetTokenT<edm::View<TGenJet> > genjetToken_;
    edm::EDGetTokenT<edm::View<TJet> > jetToken_;
    edm::EDGetTokenT<edm::View<TMET> > metToken_;
    edm::EDGetTokenT<edm::View<TMuon> > muonToken_;
    edm::EDGetTokenT<edm::View<TPhoton> > photonToken_;
    edm::EDGetTokenT<edm::View<TTau> > tauToken_;

    TH1F* hElectron_phi_;
    TH1F* hElectron_eta_;
    TH1F* hElectron_pt_;
    TH1F* hElectron_mass_;
    TH1F* hElectron_chIso_;
    TH1F* hElectron_nhIso_;
    TH1F* hElectron_phIso_;
    TH1F* hElectron_puIso_;

    TH1F* hGenJet_phi_;
    TH1F* hGenJet_eta_;
    TH1F* hGenJet_pt_;
    TH1F* hGenJet_mass_;

    TH1F* hJet_phi_;
    TH1F* hJet_eta_;
    TH1F* hJet_pt_;
    TH1F* hJet_mass_;

    TH1F* hMet_phi_;
    TH1F* hMet_eta_;
    TH1F* hMet_pt_;
    TH1F* hMet_mass_;

    TH1F* hMuon_phi_;
    TH1F* hMuon_eta_;
    TH1F* hMuon_pt_;
    TH1F* hMuon_mass_;
    TH1F* hMuon_chIso_;
    TH1F* hMuon_nhIso_;
    TH1F* hMuon_phIso_;
    TH1F* hMuon_puIso_;

    TH1F* hPhoton_phi_;
    TH1F* hPhoton_eta_;
    TH1F* hPhoton_pt_;
    TH1F* hPhoton_mass_;
    TH1F* hPhoton_chIso_;
    TH1F* hPhoton_nhIso_;
    TH1F* hPhoton_phIso_;
    TH1F* hPhoton_puIso_;

    TH1F* hTau_phi_;
    TH1F* hTau_eta_;
    TH1F* hTau_pt_;
    TH1F* hTau_mass_;

};

CATHisAnalysis::CATHisAnalysis(const edm::ParameterSet& pset )
{
  electronToken_ = consumes<edm::View<TElectron> >(pset.getParameter<edm::InputTag>("electrons"));
  genjetToken_   = consumes<edm::View<TGenJet> >(pset.getParameter<edm::InputTag>("genjets"));
  jetToken_      = consumes<edm::View<TJet> >(pset.getParameter<edm::InputTag>("jets"));
  metToken_      = consumes<edm::View<TMET> >(pset.getParameter<edm::InputTag>("mets"));
  muonToken_     = consumes<edm::View<TMuon> >(pset.getParameter<edm::InputTag>("muons"));
  photonToken_   = consumes<edm::View<TPhoton> >(pset.getParameter<edm::InputTag>("photons"));
  tauToken_      = consumes<edm::View<TTau> >(pset.getParameter<edm::InputTag>("taus"));

  edm::Service<TFileService> fs;

  TFileDirectory dirElectron = fs->mkdir("electron", "electron");
  hElectron_phi_   = dirElectron.make<TH1F>("phi","phi",100,-4,4);
  hElectron_eta_   = dirElectron.make<TH1F>("eta","eta",400,-3,3);
  hElectron_pt_    = dirElectron.make<TH1F>("pt","pt",500,0,500);
  hElectron_mass_  = dirElectron.make<TH1F>("mass","mass",1000,0,100);
  hElectron_chIso_ = dirElectron.make<TH1F>("chIso","chIso",400,0,4);
  hElectron_nhIso_ = dirElectron.make<TH1F>("nhIso","nhIso",400,0,4);
  hElectron_phIso_ = dirElectron.make<TH1F>("phIso","phIso",400,0,4);
  hElectron_puIso_ = dirElectron.make<TH1F>("puIso","puIso",400,0,4);

  TFileDirectory dirGenJet = fs->mkdir("genjet", "genjet");
  hGenJet_phi_   = dirGenJet.make<TH1F>("phi","phi",100,-4,4);
  hGenJet_eta_   = dirGenJet.make<TH1F>("eta","eta",400,-3,3);
  hGenJet_pt_    = dirGenJet.make<TH1F>("pt","pt",500,0,500);
  hGenJet_mass_  = dirGenJet.make<TH1F>("mass","mass",1000,0,100);

  TFileDirectory dirJet = fs->mkdir("jet", "jet");
  hJet_phi_   = dirJet.make<TH1F>("phi","phi",100,-4,4);
  hJet_eta_   = dirJet.make<TH1F>("eta","eta",400,-3,3);
  hJet_pt_    = dirJet.make<TH1F>("pt","pt",500,0,500);
  hJet_mass_  = dirJet.make<TH1F>("mass","mass",1000,0,100);

  TFileDirectory dirMET = fs->mkdir("MET", "MET");
  hMet_phi_   = dirMET.make<TH1F>("phi","phi",100,-4,4);
  hMet_eta_   = dirMET.make<TH1F>("eta","eta",400,-3,3);
  hMet_pt_    = dirMET.make<TH1F>("pt","pt",500,0,500);
  hMet_mass_  = dirMET.make<TH1F>("mass","mass",1000,0,100);

  TFileDirectory dirMuon = fs->mkdir("Muon", "Muon");
  hMuon_phi_   = dirMuon.make<TH1F>("phi","phi",100,-4,4);
  hMuon_eta_   = dirMuon.make<TH1F>("eta","eta",400,-3,3);
  hMuon_pt_    = dirMuon.make<TH1F>("pt","pt",500,0,500);
  hMuon_mass_  = dirMuon.make<TH1F>("mass","mass",1000,0,100);
  hMuon_chIso_ = dirMuon.make<TH1F>("chIso","chIso",400,0,4);
  hMuon_nhIso_ = dirMuon.make<TH1F>("nhIso","nhIso",400,0,4);
  hMuon_phIso_ = dirMuon.make<TH1F>("phIso","phIso",400,0,4);
  hMuon_puIso_ = dirMuon.make<TH1F>("puIso","puIso",400,0,4);

  TFileDirectory dirPhoton = fs->mkdir("photon", "photon");
  hPhoton_phi_   = dirPhoton.make<TH1F>("phi","phi",100,-4,4);
  hPhoton_eta_   = dirPhoton.make<TH1F>("eta","eta",400,-3,3);
  hPhoton_pt_    = dirPhoton.make<TH1F>("pt","pt",500,0,500);
  hPhoton_mass_  = dirPhoton.make<TH1F>("mass","mass",1000,0,100);
  hPhoton_chIso_ = dirPhoton.make<TH1F>("chIso","chIso",400,0,4);
  hPhoton_nhIso_ = dirPhoton.make<TH1F>("nhIso","nhIso",400,0,4);
  hPhoton_phIso_ = dirPhoton.make<TH1F>("phIso","phIso",400,0,4);
  hPhoton_puIso_ = dirPhoton.make<TH1F>("puIso","puIso",400,0,4);

  TFileDirectory dirTau = fs->mkdir("tau", "tau");
  hTau_phi_   = dirTau.make<TH1F>("phi","phi",100,-4,4);
  hTau_eta_   = dirTau.make<TH1F>("eta","eta",400,-3,3);
  hTau_pt_    = dirTau.make<TH1F>("pt","pt",500,0,500);
  hTau_mass_  = dirTau.make<TH1F>("mass","mass",1000,0,100);

}

void CATHisAnalysis::analyze( const edm::Event& iEvent, const edm::EventSetup&)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  edm::Handle<edm::View<TElectron> > electronHandle;
  iEvent.getByToken(electronToken_,electronHandle);
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    const auto& electron = electronHandle->at(i);

    hElectron_phi_  ->Fill( electron.phi() );
    hElectron_eta_  ->Fill( electron.eta() );
    hElectron_pt_   ->Fill( electron.pt() );
    hElectron_mass_ ->Fill( electron.mass() );
    hElectron_chIso_->Fill( electron.chargedHadronIso() );
    hElectron_nhIso_->Fill( electron.neutralHadronIso() );
    hElectron_phIso_->Fill( electron.photonIso() );
    hElectron_puIso_->Fill( electron.puChargedHadronIso() );
  }

  edm::Handle<edm::View<TGenJet> > genjetHandle;
  iEvent.getByToken(genjetToken_,genjetHandle);
  for ( int i=0, n=genjetHandle->size(); i<n; ++i )
  {
    const auto& genjet = genjetHandle->at(i);

    hGenJet_phi_  ->Fill( genjet.phi() );
    hGenJet_eta_  ->Fill( genjet.eta() );
    hGenJet_pt_   ->Fill( genjet.pt() );
    hGenJet_mass_ ->Fill( genjet.mass() );
  }

  edm::Handle<edm::View<TJet> > jetHandle;
  iEvent.getByToken(jetToken_,jetHandle);
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    const auto& jet = jetHandle->at(i);

    hJet_phi_  ->Fill( jet.phi() );
    hJet_eta_  ->Fill( jet.eta() );
    hJet_pt_   ->Fill( jet.pt() );
    hJet_mass_ ->Fill( jet.mass() );
  }

  edm::Handle<edm::View<TMET> > metHandle;
  iEvent.getByToken(metToken_,metHandle);
  for ( int i=0, n=metHandle->size(); i<n; ++i )
  {
    const auto& met = metHandle->at(i);

    hMet_phi_  ->Fill( met.phi() );
    //hMet_eta_  ->Fill( met.eta() );
    hMet_pt_   ->Fill( met.pt() );
    //hMet_mass_ ->Fill( met.mass() );
  }

  edm::Handle<edm::View<TMuon> > muonHandle;
  iEvent.getByToken(muonToken_,muonHandle);
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const auto& muon = muonHandle->at(i);

    hMuon_phi_  ->Fill( muon.phi() );
    hMuon_eta_  ->Fill( muon.eta() );
    hMuon_pt_   ->Fill( muon.pt() );
    hMuon_mass_ ->Fill( muon.mass() );
    hMuon_chIso_->Fill( muon.chargedHadronIso() );
    hMuon_nhIso_->Fill( muon.neutralHadronIso() );
    hMuon_phIso_->Fill( muon.photonIso() );
    hMuon_puIso_->Fill( muon.puChargedHadronIso() );
  }

  edm::Handle<edm::View<TPhoton> > photonHandle;
  iEvent.getByToken(photonToken_,photonHandle);
  for ( int i=0, n=photonHandle->size(); i<n; ++i )
  {
    const auto& photon = photonHandle->at(i);

    hPhoton_phi_  ->Fill( photon.phi() );
    hPhoton_eta_  ->Fill( photon.eta() );
    hPhoton_pt_   ->Fill( photon.pt() );
    hPhoton_mass_ ->Fill( photon.mass() );
    hPhoton_chIso_->Fill( photon.chargedHadronIso() );
    hPhoton_nhIso_->Fill( photon.neutralHadronIso() );
    hPhoton_phIso_->Fill( photon.photonIso() );
    hPhoton_puIso_->Fill( photon.puChargedHadronIso() );
  }

  edm::Handle<edm::View<TTau> > tauHandle;
  iEvent.getByToken(tauToken_,tauHandle);
  for ( int i=0, n=tauHandle->size(); i<n; ++i )
  {
    const auto& tau = tauHandle->at(i);

    hTau_phi_  ->Fill( tau.phi() );
    hTau_eta_  ->Fill( tau.eta() );
    hTau_pt_   ->Fill( tau.pt() );
    hTau_mass_ ->Fill( tau.mass() );
  }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATHisAnalysis);

