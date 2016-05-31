#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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

#include "TH1D.h"
#include "TTree.h"

using namespace std;

class CATHisAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit CATHisAnalysis(const edm::ParameterSet & pset);
  ~CATHisAnalysis() {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;

  edm::EDGetTokenT<cat::MuonCollection> muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<cat::METCollection> metToken_;
  edm::EDGetTokenT<cat::PhotonCollection> photonToken_;
  //edm::EDGetTokenT<cat::TauCollection> tauToken_;

  const bool isMC_;

  typedef TH1D* H1;

  // GenJets
  H1 hGenJet_n_;
  H1 hGenJet_pt_, hGenJet_eta_, hGenJet_phi_, hGenJet_mass_;
  H1 hGenJet1_pt_, hGenJet1_eta_, hGenJet1_phi_, hGenJet1_mass_;
  H1 hGenJet2_pt_, hGenJet2_eta_, hGenJet2_phi_, hGenJet2_mass_;
  H1 hGenJet3_pt_, hGenJet3_eta_, hGenJet3_phi_, hGenJet3_mass_;
  H1 hGenJet4_pt_, hGenJet4_eta_, hGenJet4_phi_, hGenJet4_mass_;

  // Electrons
  H1 hElectron_n_;
  H1 hElectron_pt_, hElectron_eta_, hElectron_phi_, hElectron_mass_;
  H1 hElectron_chIso_, hElectron_nhIso_, hElectron_phIso_, hElectron_puIso_;

  H1 hElectron1_pt_, hElectron1_eta_, hElectron1_phi_, hElectron1_mass_;
  H1 hElectron1_chIso_, hElectron1_nhIso_, hElectron1_phIso_, hElectron1_puIso_;
  H1 hElectron2_pt_, hElectron2_eta_, hElectron2_phi_, hElectron2_mass_;
  H1 hElectron2_chIso_, hElectron2_nhIso_, hElectron2_phIso_, hElectron2_puIso_;

  // Muons
  H1 hMuon_n_;
  H1 hMuon_pt_, hMuon_eta_, hMuon_phi_, hMuon_mass_;
  H1 hMuon_chIso_, hMuon_nhIso_, hMuon_phIso_, hMuon_puIso_;

  H1 hMuon1_pt_, hMuon1_eta_, hMuon1_phi_, hMuon1_mass_;
  H1 hMuon1_chIso_, hMuon1_nhIso_, hMuon1_phIso_, hMuon1_puIso_;
  H1 hMuon2_pt_, hMuon2_eta_, hMuon2_phi_, hMuon2_mass_;
  H1 hMuon2_chIso_, hMuon2_nhIso_, hMuon2_phIso_, hMuon2_puIso_;

  // Jets
  H1 hJet_n_;
  H1 hJet_pt_, hJet_eta_, hJet_phi_, hJet_mass_;
  H1 hJet_chIso_, hJet_nhIso_, hJet_phIso_, hJet_puIso_;

  H1 hJet1_pt_, hJet1_eta_, hJet1_phi_, hJet1_mass_;
  H1 hJet1_chIso_, hJet1_nhIso_, hJet1_phIso_, hJet1_puIso_;
  H1 hJet2_pt_, hJet2_eta_, hJet2_phi_, hJet2_mass_;
  H1 hJet2_chIso_, hJet2_nhIso_, hJet2_phIso_, hJet2_puIso_;
  H1 hJet3_pt_, hJet3_eta_, hJet3_phi_, hJet3_mass_;
  H1 hJet3_chIso_, hJet3_nhIso_, hJet3_phIso_, hJet3_puIso_;
  H1 hJet4_pt_, hJet4_eta_, hJet4_phi_, hJet4_mass_;
  H1 hJet4_chIso_, hJet4_nhIso_, hJet4_phIso_, hJet4_puIso_;

  // MET
  H1 hMET_pt_, hMET_phi_;

  // Photons
  H1 hPhoton_n_;
  H1 hPhoton_pt_, hPhoton_eta_, hPhoton_phi_, hPhoton_mass_;
  H1 hPhoton1_pt_, hPhoton1_eta_, hPhoton1_phi_, hPhoton1_mass_;
  H1 hPhoton2_pt_, hPhoton2_eta_, hPhoton2_phi_, hPhoton2_mass_;

  // Taus
/*
  H1 hTau_n_;
  H1 hTau_pt_, hTau_eta_, hTau_phi_, hTau_mass_;
  H1 hTau1_pt_, hTau1_eta_, hTau1_phi_, hTau1_mass_;
  H1 hTau2_pt_, hTau2_eta_, hTau2_phi_, hTau2_mass_;
*/

};

CATHisAnalysis::CATHisAnalysis(const edm::ParameterSet& pset):
  isMC_(pset.getUntrackedParameter<bool>("isMC"))
{
  if ( isMC_ ) genJetToken_   = consumes<reco::GenJetCollection>(pset.getParameter<edm::InputTag>("genJets"));
  electronToken_ = consumes<cat::ElectronCollection>(pset.getParameter<edm::InputTag>("electrons"));
  muonToken_     = consumes<cat::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));
  photonToken_   = consumes<cat::PhotonCollection>(pset.getParameter<edm::InputTag>("photons"));
  //tauToken_      = consumes<cat::TauCollection>(pset.getParameter<edm::InputTag>("taus"));
  jetToken_      = consumes<cat::JetCollection>(pset.getParameter<edm::InputTag>("jets"));
  metToken_      = consumes<cat::METCollection>(pset.getParameter<edm::InputTag>("mets"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  TFileDirectory dirElectron = fs->mkdir("electron", "electron");
  hElectron_n_    = dirElectron.make<TH1D>("h_n","n;Multiplicity",10,0,10);
  hElectron_pt_   = dirElectron.make<TH1D>("h_pt","pt;p_{T} (GeV)",500,0,500);
  hElectron_eta_  = dirElectron.make<TH1D>("h_eta","eta;#eta",400,-3,3);
  hElectron_phi_  = dirElectron.make<TH1D>("h_phi","phi;#phi",100,-4,4);
  hElectron_mass_ = dirElectron.make<TH1D>("h_mass","mass;Mass (GeV)",100,0,10);
  hElectron_chIso_ = dirElectron.make<TH1D>("h_chIso","chIso;Isolation #Delta R<0.4",400,0,4);
  hElectron_nhIso_ = dirElectron.make<TH1D>("h_nhIso","nhIso;Isolation #Delta R<0.4",400,0,4);
  hElectron_phIso_ = dirElectron.make<TH1D>("h_phIso","phIso;Isolation #Delta R<0.4",400,0,4);
  hElectron_puIso_ = dirElectron.make<TH1D>("h_puIso","puIso;Isolation #Delta R<0.4",400,0,4);

  hElectron1_pt_   = dirElectron.make<TH1D>("h1_pt","pt",500,0,500);
  hElectron1_eta_  = dirElectron.make<TH1D>("h1_eta","eta",400,-3,3);
  hElectron1_phi_  = dirElectron.make<TH1D>("h1_phi","phi",100,-4,4);
  hElectron1_mass_ = dirElectron.make<TH1D>("h1_mass","mass",100,0,10);
  hElectron1_chIso_ = dirElectron.make<TH1D>("h1_chIso","chIso",400,0,4);
  hElectron1_nhIso_ = dirElectron.make<TH1D>("h1_nhIso","nhIso",400,0,4);
  hElectron1_phIso_ = dirElectron.make<TH1D>("h1_phIso","phIso",400,0,4);
  hElectron1_puIso_ = dirElectron.make<TH1D>("h1_puIso","puIso",400,0,4);
  hElectron2_pt_   = dirElectron.make<TH1D>("h2_pt","pt",500,0,500);
  hElectron2_eta_  = dirElectron.make<TH1D>("h2_eta","eta",400,-3,3);
  hElectron2_phi_  = dirElectron.make<TH1D>("h2_phi","phi",100,-4,4);
  hElectron2_mass_ = dirElectron.make<TH1D>("h2_mass","mass",100,0,10);
  hElectron2_chIso_ = dirElectron.make<TH1D>("h2_chIso","chIso",400,0,4);
  hElectron2_nhIso_ = dirElectron.make<TH1D>("h2_nhIso","nhIso",400,0,4);
  hElectron2_phIso_ = dirElectron.make<TH1D>("h2_phIso","phIso",400,0,4);
  hElectron2_puIso_ = dirElectron.make<TH1D>("h2_puIso","puIso",400,0,4);

  TFileDirectory dirMuon = fs->mkdir("electron", "electron");
  hMuon_n_    = dirMuon.make<TH1D>("h_n","n;Multiplicity",10,0,10);
  hMuon_pt_   = dirMuon.make<TH1D>("h_pt","pt",500,0,500);
  hMuon_eta_  = dirMuon.make<TH1D>("h_eta","eta",400,-3,3);
  hMuon_phi_  = dirMuon.make<TH1D>("h_phi","phi",100,-4,4);
  hMuon_mass_ = dirMuon.make<TH1D>("h_mass","mass",100,0,10);
  hMuon_chIso_ = dirMuon.make<TH1D>("h_chIso","chIso",400,0,4);
  hMuon_nhIso_ = dirMuon.make<TH1D>("h_nhIso","nhIso",400,0,4);
  hMuon_phIso_ = dirMuon.make<TH1D>("h_phIso","phIso",400,0,4);
  hMuon_puIso_ = dirMuon.make<TH1D>("h_puIso","puIso",400,0,4);

  hMuon1_pt_   = dirMuon.make<TH1D>("h1_pt","pt",500,0,500);
  hMuon1_eta_  = dirMuon.make<TH1D>("h1_eta","eta",400,-3,3);
  hMuon1_phi_  = dirMuon.make<TH1D>("h1_phi","phi",100,-4,4);
  hMuon1_mass_ = dirMuon.make<TH1D>("h1_mass","mass",100,0,10);
  hMuon1_chIso_ = dirMuon.make<TH1D>("h1_chIso","chIso",400,0,4);
  hMuon1_nhIso_ = dirMuon.make<TH1D>("h1_nhIso","nhIso",400,0,4);
  hMuon1_phIso_ = dirMuon.make<TH1D>("h1_phIso","phIso",400,0,4);
  hMuon1_puIso_ = dirMuon.make<TH1D>("h1_puIso","puIso",400,0,4);
  hMuon2_pt_   = dirMuon.make<TH1D>("h2_pt","pt",500,0,500);
  hMuon2_eta_  = dirMuon.make<TH1D>("h2_eta","eta",400,-3,3);
  hMuon2_phi_  = dirMuon.make<TH1D>("h2_phi","phi",100,-4,4);
  hMuon2_mass_ = dirMuon.make<TH1D>("h2_mass","mass",100,0,10);
  hMuon2_chIso_ = dirMuon.make<TH1D>("h2_chIso","chIso",400,0,4);
  hMuon2_nhIso_ = dirMuon.make<TH1D>("h2_nhIso","nhIso",400,0,4);
  hMuon2_phIso_ = dirMuon.make<TH1D>("h2_phIso","phIso",400,0,4);
  hMuon2_puIso_ = dirMuon.make<TH1D>("h2_puIso","puIso",400,0,4);

/*
  TFileDirectory dirTau = fs->mkdir("tau", "tau");
  hTau_n_    = dirTau.make<TH1D>("h_n","n;Multiplicity",10,0,10);
  hTau_pt_   = dirTau.make<TH1D>("h_pt","pt",500,0,500);
  hTau_eta_  = dirTau.make<TH1D>("h_eta","eta",400,-3,3);
  hTau_phi_  = dirTau.make<TH1D>("h_phi","phi",100,-4,4);
  hTau_mass_ = dirTau.make<TH1D>("h_mass","mass",100,0,10);

  hTau1_pt_   = dirTau.make<TH1D>("h1_pt","pt",500,0,500);
  hTau1_eta_  = dirTau.make<TH1D>("h1_eta","eta",400,-3,3);
  hTau1_phi_  = dirTau.make<TH1D>("h1_phi","phi",100,-4,4);
  hTau1_mass_ = dirTau.make<TH1D>("h1_mass","mass",100,0,10);
  hTau2_pt_   = dirTau.make<TH1D>("h2_pt","pt",500,0,500);
  hTau2_eta_  = dirTau.make<TH1D>("h2_eta","eta",400,-3,3);
  hTau2_phi_  = dirTau.make<TH1D>("h2_phi","phi",100,-4,4);
  hTau2_mass_ = dirTau.make<TH1D>("h2_mass","mass",100,0,10);
*/

  TFileDirectory dirPhoton = fs->mkdir("photon", "photon");
  hPhoton_n_    = dirPhoton.make<TH1D>("h_n","n;Multiplicity",10,0,10);
  hPhoton_pt_   = dirPhoton.make<TH1D>("h_pt","pt",500,0,500);
  hPhoton_eta_  = dirPhoton.make<TH1D>("h_eta","eta",400,-3,3);
  hPhoton_phi_  = dirPhoton.make<TH1D>("h_phi","phi",100,-4,4);
  hPhoton_mass_ = dirPhoton.make<TH1D>("h_mass","mass",100,0,10);

  hPhoton1_pt_   = dirPhoton.make<TH1D>("h1_pt","pt",500,0,500);
  hPhoton1_eta_  = dirPhoton.make<TH1D>("h1_eta","eta",400,-3,3);
  hPhoton1_phi_  = dirPhoton.make<TH1D>("h1_phi","phi",100,-4,4);
  hPhoton1_mass_ = dirPhoton.make<TH1D>("h1_mass","mass",100,0,10);
  hPhoton2_pt_   = dirPhoton.make<TH1D>("h2_pt","pt",500,0,500);
  hPhoton2_eta_  = dirPhoton.make<TH1D>("h2_eta","eta",400,-3,3);
  hPhoton2_phi_  = dirPhoton.make<TH1D>("h2_phi","phi",100,-4,4);
  hPhoton2_mass_ = dirPhoton.make<TH1D>("h2_mass","mass",100,0,10);

  TFileDirectory dirJet = fs->mkdir("jet", "jet");
  hJet_n_    = dirJet.make<TH1D>("h_n","n;Multiplicity",10,0,10);
  hJet_pt_    = dirJet.make<TH1D>("h_pt","pt",500,0,500);
  hJet_eta_   = dirJet.make<TH1D>("h_eta","eta",400,-3,3);
  hJet_phi_   = dirJet.make<TH1D>("h_phi","phi",100,-4,4);
  hJet_mass_  = dirJet.make<TH1D>("h_mass","mass",1000,0,100);

  hJet1_pt_    = dirJet.make<TH1D>("h1_pt","pt",500,0,500);
  hJet1_eta_   = dirJet.make<TH1D>("h1_eta","eta",400,-3,3);
  hJet1_phi_   = dirJet.make<TH1D>("h1_phi","phi",100,-4,4);
  hJet1_mass_  = dirJet.make<TH1D>("h1_mass","mass",1000,0,100);
  hJet2_pt_    = dirJet.make<TH1D>("h2_pt","pt",500,0,500);
  hJet2_eta_   = dirJet.make<TH1D>("h2_eta","eta",400,-3,3);
  hJet2_phi_   = dirJet.make<TH1D>("h2_phi","phi",100,-4,4);
  hJet2_mass_  = dirJet.make<TH1D>("h2_mass","mass",1000,0,100);
  hJet3_pt_    = dirJet.make<TH1D>("h3_pt","pt",500,0,500);
  hJet3_eta_   = dirJet.make<TH1D>("h3_eta","eta",400,-3,3);
  hJet3_phi_   = dirJet.make<TH1D>("h3_phi","phi",100,-4,4);
  hJet3_mass_  = dirJet.make<TH1D>("h3_mass","mass",1000,0,100);
  hJet4_pt_    = dirJet.make<TH1D>("h4_pt","pt",500,0,500);
  hJet4_eta_   = dirJet.make<TH1D>("h4_eta","eta",400,-3,3);
  hJet4_phi_   = dirJet.make<TH1D>("h4_phi","phi",100,-4,4);
  hJet4_mass_  = dirJet.make<TH1D>("h4_mass","mass",1000,0,100);

  if ( isMC_ ) {
    TFileDirectory dirGenJet = fs->mkdir("genJet", "genJet");
    hGenJet_n_    = dirGenJet.make<TH1D>("h_n","n;Multiplicity",10,0,10);
    hGenJet_pt_   = dirGenJet.make<TH1D>("h_pt","pt",500,0,500);
    hGenJet_eta_  = dirGenJet.make<TH1D>("h_eta","eta",400,-3,3);
    hGenJet_phi_  = dirGenJet.make<TH1D>("h_phi","phi",100,-4,4);
    hGenJet_mass_ = dirGenJet.make<TH1D>("h_mass","mass",1000,0,100);

    hGenJet1_pt_    = dirGenJet.make<TH1D>("h1_pt","pt",500,0,500);
    hGenJet1_eta_   = dirGenJet.make<TH1D>("h1_eta","eta",400,-3,3);
    hGenJet1_phi_   = dirGenJet.make<TH1D>("h1_phi","phi",100,-4,4);
    hGenJet1_mass_  = dirGenJet.make<TH1D>("h1_mass","mass",1000,0,100);
    hGenJet2_pt_    = dirGenJet.make<TH1D>("h2_pt","pt",500,0,500);
    hGenJet2_eta_   = dirGenJet.make<TH1D>("h2_eta","eta",400,-3,3);
    hGenJet2_phi_   = dirGenJet.make<TH1D>("h2_phi","phi",100,-4,4);
    hGenJet2_mass_  = dirGenJet.make<TH1D>("h2_mass","mass",1000,0,100);
    hGenJet3_pt_    = dirGenJet.make<TH1D>("h3_pt","pt",500,0,500);
    hGenJet3_eta_   = dirGenJet.make<TH1D>("h3_eta","eta",400,-3,3);
    hGenJet3_phi_   = dirGenJet.make<TH1D>("h3_phi","phi",100,-4,4);
    hGenJet3_mass_  = dirGenJet.make<TH1D>("h3_mass","mass",1000,0,100);
    hGenJet4_pt_    = dirGenJet.make<TH1D>("h4_pt","pt",500,0,500);
    hGenJet4_eta_   = dirGenJet.make<TH1D>("h4_eta","eta",400,-3,3);
    hGenJet4_phi_   = dirGenJet.make<TH1D>("h4_phi","phi",100,-4,4);
    hGenJet4_mass_  = dirGenJet.make<TH1D>("h4_mass","mass",1000,0,100);
  }

  TFileDirectory dirMET = fs->mkdir("MET", "MET");
  hMET_pt_    = dirMET.make<TH1D>("h_pt","pt",500,0,500);
  hMET_phi_   = dirMET.make<TH1D>("h_phi","phi",100,-4,4);

}

void CATHisAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup&)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  edm::Handle<cat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken_,electronHandle);
  const int nElectron = electronHandle->size();
  hElectron_n_->Fill(nElectron);
  for ( int i=0; i<nElectron; ++i ) {
    const auto& electron = electronHandle->at(i);

    hElectron_pt_   ->Fill( electron.pt() );
    hElectron_eta_  ->Fill( electron.eta() );
    hElectron_phi_  ->Fill( electron.phi() );
    hElectron_mass_ ->Fill( electron.mass() );
    hElectron_chIso_->Fill( electron.chargedHadronIso() );
    hElectron_nhIso_->Fill( electron.neutralHadronIso() );
    hElectron_phIso_->Fill( electron.photonIso() );
    hElectron_puIso_->Fill( electron.puChargedHadronIso() );
  }
  if ( nElectron >= 1 ) {
    const auto& electron = electronHandle->at(0);
    hElectron1_pt_   ->Fill( electron.pt() );
    hElectron1_eta_  ->Fill( electron.eta() );
    hElectron1_phi_  ->Fill( electron.phi() );
    hElectron1_mass_ ->Fill( electron.mass() );
    hElectron1_chIso_->Fill( electron.chargedHadronIso() );
    hElectron1_nhIso_->Fill( electron.neutralHadronIso() );
    hElectron1_phIso_->Fill( electron.photonIso() );
    hElectron1_puIso_->Fill( electron.puChargedHadronIso() );
  }
  if ( nElectron >= 2 ) {
    const auto& electron = electronHandle->at(1);
    hElectron2_pt_   ->Fill( electron.pt() );
    hElectron2_eta_  ->Fill( electron.eta() );
    hElectron2_phi_  ->Fill( electron.phi() );
    hElectron2_mass_ ->Fill( electron.mass() );
    hElectron2_chIso_->Fill( electron.chargedHadronIso() );
    hElectron2_nhIso_->Fill( electron.neutralHadronIso() );
    hElectron2_phIso_->Fill( electron.photonIso() );
    hElectron2_puIso_->Fill( electron.puChargedHadronIso() );
  }

  edm::Handle<cat::MuonCollection> muonHandle;
  iEvent.getByToken(muonToken_,muonHandle);
  const int nMuon = muonHandle->size();
  hMuon_n_->Fill(nMuon);
  for ( int i=0; i<nMuon; ++i ) {
    const auto& muon = muonHandle->at(i);

    hMuon_pt_   ->Fill( muon.pt() );
    hMuon_eta_  ->Fill( muon.eta() );
    hMuon_phi_  ->Fill( muon.phi() );
    hMuon_mass_ ->Fill( muon.mass() );
    hMuon_chIso_->Fill( muon.chargedHadronIso() );
    hMuon_nhIso_->Fill( muon.neutralHadronIso() );
    hMuon_phIso_->Fill( muon.photonIso() );
    hMuon_puIso_->Fill( muon.puChargedHadronIso() );
  }
  if ( nMuon >= 1 ) {
    const auto& muon = muonHandle->at(0);
    hMuon1_pt_   ->Fill( muon.pt() );
    hMuon1_eta_  ->Fill( muon.eta() );
    hMuon1_phi_  ->Fill( muon.phi() );
    hMuon1_mass_ ->Fill( muon.mass() );
    hMuon1_chIso_->Fill( muon.chargedHadronIso() );
    hMuon1_nhIso_->Fill( muon.neutralHadronIso() );
    hMuon1_phIso_->Fill( muon.photonIso() );
    hMuon1_puIso_->Fill( muon.puChargedHadronIso() );
  }
  if ( nMuon >= 2 ) {
    const auto& muon = muonHandle->at(1);
    hMuon2_pt_   ->Fill( muon.pt() );
    hMuon2_eta_  ->Fill( muon.eta() );
    hMuon2_phi_  ->Fill( muon.phi() );
    hMuon2_mass_ ->Fill( muon.mass() );
    hMuon2_chIso_->Fill( muon.chargedHadronIso() );
    hMuon2_nhIso_->Fill( muon.neutralHadronIso() );
    hMuon2_phIso_->Fill( muon.photonIso() );
    hMuon2_puIso_->Fill( muon.puChargedHadronIso() );
  }

  edm::Handle<cat::PhotonCollection> photonHandle;
  iEvent.getByToken(photonToken_,photonHandle);
  const int nPhoton = photonHandle->size();
  hPhoton_n_->Fill(nPhoton);
  for ( int i=0; i<nPhoton; ++i ) {
    const auto& photon = photonHandle->at(i);

    hPhoton_pt_   ->Fill( photon.pt() );
    hPhoton_eta_  ->Fill( photon.eta() );
    hPhoton_phi_  ->Fill( photon.phi() );
    hPhoton_mass_ ->Fill( photon.mass() );
  }
  if ( nPhoton >= 1 ) {
    const auto& photon = photonHandle->at(0);
    hPhoton1_pt_   ->Fill( photon.pt() );
    hPhoton1_eta_  ->Fill( photon.eta() );
    hPhoton1_phi_  ->Fill( photon.phi() );
    hPhoton1_mass_ ->Fill( photon.mass() );
  }
  if ( nPhoton >= 2 ) {
    const auto& photon = photonHandle->at(1);
    hPhoton2_pt_   ->Fill( photon.pt() );
    hPhoton2_eta_  ->Fill( photon.eta() );
    hPhoton2_phi_  ->Fill( photon.phi() );
    hPhoton2_mass_ ->Fill( photon.mass() );
  }

/*
  edm::Handle<cat::TauCollection> tauHandle;
  iEvent.getByToken(tauToken_,tauHandle);
  const int nTau = tauHandle->size();
  hTau_n_->Fill(nTau);
  for ( int i=0; i<nTau; ++i ) {
    const auto& tau = tauHandle->at(i);

    hTau_pt_   ->Fill( tau.pt() );
    hTau_eta_  ->Fill( tau.eta() );
    hTau_phi_  ->Fill( tau.phi() );
    hTau_mass_ ->Fill( tau.mass() );
  }
  if ( nTau >= 1 ) {
    const auto& tau = tauHandle->at(0);
    hTau1_pt_   ->Fill( tau.pt() );
    hTau1_eta_  ->Fill( tau.eta() );
    hTau1_phi_  ->Fill( tau.phi() );
    hTau1_mass_ ->Fill( tau.mass() );
  }
  if ( nTau >= 2 ) {
    const auto& tau = tauHandle->at(1);
    hTau2_pt_   ->Fill( tau.pt() );
    hTau2_eta_  ->Fill( tau.eta() );
    hTau2_phi_  ->Fill( tau.phi() );
    hTau2_mass_ ->Fill( tau.mass() );
  }
*/

  edm::Handle<cat::JetCollection> jetHandle;
  iEvent.getByToken(jetToken_,jetHandle);
  const int nJet = jetHandle->size();
  hJet_n_->Fill(nJet);
  for ( int i=0; i<nJet; ++i ) {
    const auto& jet = jetHandle->at(i);

    hJet_pt_   ->Fill( jet.pt() );
    hJet_eta_  ->Fill( jet.eta() );
    hJet_phi_  ->Fill( jet.phi() );
    hJet_mass_ ->Fill( jet.mass() );
  }
  if ( nJet >= 1 ) {
    const auto& jet = jetHandle->at(0);
    hJet1_pt_   ->Fill( jet.pt() );
    hJet1_eta_  ->Fill( jet.eta() );
    hJet1_phi_  ->Fill( jet.phi() );
    hJet1_mass_ ->Fill( jet.mass() );
  }
  if ( nJet >= 2 ) {
    const auto& jet = jetHandle->at(1);
    hJet2_pt_   ->Fill( jet.pt() );
    hJet2_eta_  ->Fill( jet.eta() );
    hJet2_phi_  ->Fill( jet.phi() );
    hJet2_mass_ ->Fill( jet.mass() );
  }
  if ( nJet >= 3 ) {
    const auto& jet = jetHandle->at(2);
    hJet3_pt_   ->Fill( jet.pt() );
    hJet3_eta_  ->Fill( jet.eta() );
    hJet3_phi_  ->Fill( jet.phi() );
    hJet3_mass_ ->Fill( jet.mass() );
  }
  if ( nJet >= 4 ) {
    const auto& jet = jetHandle->at(3);
    hJet3_pt_   ->Fill( jet.pt() );
    hJet3_eta_  ->Fill( jet.eta() );
    hJet3_phi_  ->Fill( jet.phi() );
    hJet3_mass_ ->Fill( jet.mass() );
  }

  if ( isMC_ ) {
    edm::Handle<reco::GenJetCollection> genJetHandle;
    iEvent.getByToken(genJetToken_,genJetHandle);
    const int nGenJet = genJetHandle->size();
    hGenJet_n_->Fill(nGenJet);
    for ( int i=0; i<nGenJet; ++i ) {
      const auto& genJet = genJetHandle->at(i);

      hGenJet_pt_   ->Fill( genJet.pt() );
      hGenJet_eta_  ->Fill( genJet.eta() );
      hGenJet_phi_  ->Fill( genJet.phi() );
      hGenJet_mass_ ->Fill( genJet.mass() );
    }
    if ( nGenJet >= 1 ) {
      const auto& genJet = genJetHandle->at(0);
      hGenJet1_pt_   ->Fill( genJet.pt() );
      hGenJet1_eta_  ->Fill( genJet.eta() );
      hGenJet1_phi_  ->Fill( genJet.phi() );
      hGenJet1_mass_ ->Fill( genJet.mass() );
    }
    if ( nGenJet >= 2 ) {
      const auto& genJet = genJetHandle->at(1);
      hGenJet2_pt_   ->Fill( genJet.pt() );
      hGenJet2_eta_  ->Fill( genJet.eta() );
      hGenJet2_phi_  ->Fill( genJet.phi() );
      hGenJet2_mass_ ->Fill( genJet.mass() );
    }
    if ( nGenJet >= 3 ) {
      const auto& genJet = genJetHandle->at(2);
      hGenJet3_pt_   ->Fill( genJet.pt() );
      hGenJet3_eta_  ->Fill( genJet.eta() );
      hGenJet3_phi_  ->Fill( genJet.phi() );
      hGenJet3_mass_ ->Fill( genJet.mass() );
    }
    if ( nGenJet >= 4 ) {
      const auto& genJet = genJetHandle->at(3);
      hGenJet3_pt_   ->Fill( genJet.pt() );
      hGenJet3_eta_  ->Fill( genJet.eta() );
      hGenJet3_phi_  ->Fill( genJet.phi() );
      hGenJet3_mass_ ->Fill( genJet.mass() );
    }
  }

  edm::Handle<cat::METCollection> metHandle;
  iEvent.getByToken(metToken_,metHandle);
  if ( !metHandle->empty() ) {
    const auto& met = metHandle->at(0);
    hMET_pt_   ->Fill( met.pt() );
    hMET_phi_  ->Fill( met.phi() );
  }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATHisAnalysis);

