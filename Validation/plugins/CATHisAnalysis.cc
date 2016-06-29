#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/View.h"

#include "CATTools/DataFormats/interface/GenWeights.h"
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
  edm::EDGetTokenT<int> nPVToken_;
  edm::EDGetTokenT<float> puWeightToken_, genWeightToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;

  edm::EDGetTokenT<cat::MuonCollection> muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<cat::METCollection> metToken_;
  edm::EDGetTokenT<cat::PhotonCollection> photonToken_;
  //edm::EDGetTokenT<cat::TauCollection> tauToken_;

  const bool isMC_;

  typedef TH1D* H1;

  H1 hNPV_, hNPVRaw_;

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
  H1 hElectron_ids_, hElectron_enUp_, hElectron_enDn_;

  H1 hElectron1_pt_, hElectron1_eta_, hElectron1_phi_, hElectron1_mass_;
  H1 hElectron1_chIso_, hElectron1_nhIso_, hElectron1_phIso_, hElectron1_puIso_;
  H1 hElectron1_ids_, hElectron1_enUp_, hElectron1_enDn_;
  H1 hElectron2_pt_, hElectron2_eta_, hElectron2_phi_, hElectron2_mass_;
  H1 hElectron2_chIso_, hElectron2_nhIso_, hElectron2_phIso_, hElectron2_puIso_;
  H1 hElectron2_ids_, hElectron2_enUp_, hElectron2_enDn_;

  // Muons
  H1 hMuon_n_;
  H1 hMuon_pt_, hMuon_eta_, hMuon_phi_, hMuon_mass_;
  H1 hMuon_chIso_, hMuon_nhIso_, hMuon_phIso_, hMuon_puIso_;
  H1 hMuon_ids_, hMuon_enUp_, hMuon_enDn_;

  H1 hMuon1_pt_, hMuon1_eta_, hMuon1_phi_, hMuon1_mass_;
  H1 hMuon1_chIso_, hMuon1_nhIso_, hMuon1_phIso_, hMuon1_puIso_;
  H1 hMuon1_ids_, hMuon1_enUp_, hMuon1_enDn_;
  H1 hMuon2_pt_, hMuon2_eta_, hMuon2_phi_, hMuon2_mass_;
  H1 hMuon2_chIso_, hMuon2_nhIso_, hMuon2_phIso_, hMuon2_puIso_;
  H1 hMuon2_ids_, hMuon2_enUp_, hMuon2_enDn_;

  // Jets
  H1 hJet_n_;
  H1 hJet_pt_, hJet_eta_, hJet_phi_, hJet_mass_;
  H1 hJet_chIso_, hJet_nhIso_, hJet_phIso_, hJet_puIso_;
  H1 hJet_jesUp_, hJet_jesDn_, hJet_jer_, hJet_jerUp_, hJet_jerDn_, hJet_CSV_;

  H1 hJet1_pt_, hJet1_eta_, hJet1_phi_, hJet1_mass_;
  H1 hJet1_chIso_, hJet1_nhIso_, hJet1_phIso_, hJet1_puIso_;
  H1 hJet1_jesUp_, hJet1_jesDn_, hJet1_jer_, hJet1_jerUp_, hJet1_jerDn_, hJet1_CSV_;
  H1 hJet2_pt_, hJet2_eta_, hJet2_phi_, hJet2_mass_;
  H1 hJet2_chIso_, hJet2_nhIso_, hJet2_phIso_, hJet2_puIso_;
  H1 hJet2_jesUp_, hJet2_jesDn_, hJet2_jer_, hJet2_jerUp_, hJet2_jerDn_, hJet2_CSV_;
  H1 hJet3_pt_, hJet3_eta_, hJet3_phi_, hJet3_mass_;
  H1 hJet3_chIso_, hJet3_nhIso_, hJet3_phIso_, hJet3_puIso_;
  H1 hJet3_jesUp_, hJet3_jesDn_, hJet3_jer_, hJet3_jerUp_, hJet3_jerDn_, hJet3_CSV_;
  H1 hJet4_pt_, hJet4_eta_, hJet4_phi_, hJet4_mass_;
  H1 hJet4_chIso_, hJet4_nhIso_, hJet4_phIso_, hJet4_puIso_;
  H1 hJet4_jesUp_, hJet4_jesDn_, hJet4_jer_, hJet4_jerUp_, hJet4_jerDn_, hJet4_CSV_;

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
  if ( isMC_ ) {
    genWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("genWeight"));
    puWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("puWeight"));
    genJetToken_   = consumes<reco::GenJetCollection>(pset.getParameter<edm::InputTag>("genJets"));
  }
  nPVToken_ = consumes<int>(pset.getParameter<edm::InputTag>("nPV"));
  electronToken_ = consumes<cat::ElectronCollection>(pset.getParameter<edm::InputTag>("electrons"));
  muonToken_     = consumes<cat::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));
  photonToken_   = consumes<cat::PhotonCollection>(pset.getParameter<edm::InputTag>("photons"));
  //tauToken_      = consumes<cat::TauCollection>(pset.getParameter<edm::InputTag>("taus"));
  jetToken_      = consumes<cat::JetCollection>(pset.getParameter<edm::InputTag>("jets"));
  metToken_      = consumes<cat::METCollection>(pset.getParameter<edm::InputTag>("mets"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  hNPV_ = fs->make<TH1D>("hNPV", "nPV;Vertex multiplicity;Events", 50, 0, 50);
  hNPVRaw_ = fs->make<TH1D>("hNPVRaw", "nPV no weight;Vertex multiplicity;Events", 50, 0, 50);

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
  hElectron_ids_  = dirElectron.make<TH1D>("h_id","id;ids",8,1,9);
  hElectron_ids_->GetXaxis()->SetBinLabel(1, "Veto");
  hElectron_ids_->GetXaxis()->SetBinLabel(2, "MediumMVA");
  hElectron_ids_->GetXaxis()->SetBinLabel(3, "TightMVA");
  hElectron_ids_->GetXaxis()->SetBinLabel(4, "CutBased-loose");
  hElectron_ids_->GetXaxis()->SetBinLabel(5, "CutBased-medium");
  hElectron_ids_->GetXaxis()->SetBinLabel(6, "CutBased-tight");
  hElectron_ids_->GetXaxis()->SetBinLabel(7, "MVA-WP80");
  hElectron_ids_->GetXaxis()->SetBinLabel(8, "MVA-WP90");
  hElectron_enUp_ = dirElectron.make<TH1D>("h_enUp","enUp;Energy scale up",100,0,3);
  hElectron_enDn_ = dirElectron.make<TH1D>("h_enDn","enDn;Energy scale down",100,0,3);

  hElectron1_pt_   = dirElectron.make<TH1D>("h1_pt","pt",500,0,500);
  hElectron1_eta_  = dirElectron.make<TH1D>("h1_eta","eta",400,-3,3);
  hElectron1_phi_  = dirElectron.make<TH1D>("h1_phi","phi",100,-4,4);
  hElectron1_mass_ = dirElectron.make<TH1D>("h1_mass","mass",100,0,10);
  hElectron1_chIso_ = dirElectron.make<TH1D>("h1_chIso","chIso",400,0,4);
  hElectron1_nhIso_ = dirElectron.make<TH1D>("h1_nhIso","nhIso",400,0,4);
  hElectron1_phIso_ = dirElectron.make<TH1D>("h1_phIso","phIso",400,0,4);
  hElectron1_puIso_ = dirElectron.make<TH1D>("h1_puIso","puIso",400,0,4);
  hElectron1_ids_  = dirElectron.make<TH1D>("h1_id","id;ids",8,1,9);
  hElectron1_ids_->GetXaxis()->SetBinLabel(1, "Veto");
  hElectron1_ids_->GetXaxis()->SetBinLabel(2, "MediumMVA");
  hElectron1_ids_->GetXaxis()->SetBinLabel(3, "TightMVA");
  hElectron1_ids_->GetXaxis()->SetBinLabel(4, "CutBased-loose");
  hElectron1_ids_->GetXaxis()->SetBinLabel(5, "CutBased-medium");
  hElectron1_ids_->GetXaxis()->SetBinLabel(6, "CutBased-tight");
  hElectron1_ids_->GetXaxis()->SetBinLabel(7, "MVA-WP80");
  hElectron1_ids_->GetXaxis()->SetBinLabel(8, "MVA-WP90");
  hElectron1_enUp_ = dirElectron.make<TH1D>("h1_enUp","enUp;Energy scale up",100,0,3);
  hElectron1_enDn_ = dirElectron.make<TH1D>("h1_enDn","enDn;Energy scale down",100,0,3);

  hElectron2_pt_   = dirElectron.make<TH1D>("h2_pt","pt",500,0,500);
  hElectron2_eta_  = dirElectron.make<TH1D>("h2_eta","eta",400,-3,3);
  hElectron2_phi_  = dirElectron.make<TH1D>("h2_phi","phi",100,-4,4);
  hElectron2_mass_ = dirElectron.make<TH1D>("h2_mass","mass",100,0,10);
  hElectron2_chIso_ = dirElectron.make<TH1D>("h2_chIso","chIso",400,0,4);
  hElectron2_nhIso_ = dirElectron.make<TH1D>("h2_nhIso","nhIso",400,0,4);
  hElectron2_phIso_ = dirElectron.make<TH1D>("h2_phIso","phIso",400,0,4);
  hElectron2_puIso_ = dirElectron.make<TH1D>("h2_puIso","puIso",400,0,4);
  hElectron2_ids_  = dirElectron.make<TH1D>("h2_id","id;ids",8,1,9);
  hElectron2_ids_->GetXaxis()->SetBinLabel(1, "Veto");
  hElectron2_ids_->GetXaxis()->SetBinLabel(2, "MediumMVA");
  hElectron2_ids_->GetXaxis()->SetBinLabel(3, "TightMVA");
  hElectron2_ids_->GetXaxis()->SetBinLabel(4, "CutBased-loose");
  hElectron2_ids_->GetXaxis()->SetBinLabel(5, "CutBased-medium");
  hElectron2_ids_->GetXaxis()->SetBinLabel(6, "CutBased-tight");
  hElectron2_ids_->GetXaxis()->SetBinLabel(7, "MVA-WP80");
  hElectron2_ids_->GetXaxis()->SetBinLabel(8, "MVA-WP90");
  hElectron2_enUp_ = dirElectron.make<TH1D>("h2_enUp","enUp;Energy scale up",100,0,3);
  hElectron2_enDn_ = dirElectron.make<TH1D>("h2_enDn","enDn;Energy scale down",100,0,3);

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
  hMuon_ids_  = dirMuon.make<TH1D>("h_id","id;ids",2,1,3);
  hMuon_ids_->GetXaxis()->SetBinLabel(1, "Loose");
  hMuon_ids_->GetXaxis()->SetBinLabel(2, "Tight");
  hMuon_enUp_ = dirMuon.make<TH1D>("h_enUp","enUp;Energy scale up",100,0,3);
  hMuon_enDn_ = dirMuon.make<TH1D>("h_enDn","enDn;Energy scale down",100,0,3);

  hMuon1_pt_   = dirMuon.make<TH1D>("h1_pt","pt",500,0,500);
  hMuon1_eta_  = dirMuon.make<TH1D>("h1_eta","eta",400,-3,3);
  hMuon1_phi_  = dirMuon.make<TH1D>("h1_phi","phi",100,-4,4);
  hMuon1_mass_ = dirMuon.make<TH1D>("h1_mass","mass",100,0,10);
  hMuon1_chIso_ = dirMuon.make<TH1D>("h1_chIso","chIso",400,0,4);
  hMuon1_nhIso_ = dirMuon.make<TH1D>("h1_nhIso","nhIso",400,0,4);
  hMuon1_phIso_ = dirMuon.make<TH1D>("h1_phIso","phIso",400,0,4);
  hMuon1_puIso_ = dirMuon.make<TH1D>("h1_puIso","puIso",400,0,4);
  hMuon1_ids_  = dirMuon.make<TH1D>("h1_id","id;ids",2,1,3);
  hMuon1_ids_->GetXaxis()->SetBinLabel(1, "Loose");
  hMuon1_ids_->GetXaxis()->SetBinLabel(2, "Tight");
  hMuon1_enUp_ = dirMuon.make<TH1D>("h1_enUp","enUp;Energy scale up",100,0,3);
  hMuon1_enDn_ = dirMuon.make<TH1D>("h1_enDn","enDn;Energy scale down",100,0,3);
  hMuon2_pt_   = dirMuon.make<TH1D>("h2_pt","pt",500,0,500);
  hMuon2_eta_  = dirMuon.make<TH1D>("h2_eta","eta",400,-3,3);
  hMuon2_phi_  = dirMuon.make<TH1D>("h2_phi","phi",100,-4,4);
  hMuon2_mass_ = dirMuon.make<TH1D>("h2_mass","mass",100,0,10);
  hMuon2_chIso_ = dirMuon.make<TH1D>("h2_chIso","chIso",400,0,4);
  hMuon2_nhIso_ = dirMuon.make<TH1D>("h2_nhIso","nhIso",400,0,4);
  hMuon2_phIso_ = dirMuon.make<TH1D>("h2_phIso","phIso",400,0,4);
  hMuon2_puIso_ = dirMuon.make<TH1D>("h2_puIso","puIso",400,0,4);
  hMuon2_ids_  = dirMuon.make<TH1D>("h2_id","id;ids",2,1,3);
  hMuon2_ids_->GetXaxis()->SetBinLabel(1, "Loose");
  hMuon2_ids_->GetXaxis()->SetBinLabel(2, "Tight");
  hMuon2_enUp_ = dirMuon.make<TH1D>("h2_enUp","enUp;Energy scale up",100,0,3);
  hMuon2_enDn_ = dirMuon.make<TH1D>("h2_enDn","enDn;Energy scale down",100,0,3);

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
  hJet_jesUp_ = dirJet.make<TH1D>("h_jesUp","jesUp",100,0,3);
  hJet_jesDn_ = dirJet.make<TH1D>("h_jesDn","jesDn",100,0,3);
  hJet_jer_   = dirJet.make<TH1D>("h_jer"  ,"jer"  ,100,0,3);
  hJet_jerUp_ = dirJet.make<TH1D>("h_jerUp","jerUp",100,0,3);
  hJet_jerDn_ = dirJet.make<TH1D>("h_jerDn","jerDn",100,0,3);
  hJet_CSV_ = dirJet.make<TH1D>("h_CSV","CSV",100,-1,1);

  hJet1_pt_    = dirJet.make<TH1D>("h1_pt","pt",500,0,500);
  hJet1_eta_   = dirJet.make<TH1D>("h1_eta","eta",400,-3,3);
  hJet1_phi_   = dirJet.make<TH1D>("h1_phi","phi",100,-4,4);
  hJet1_mass_  = dirJet.make<TH1D>("h1_mass","mass",1000,0,100);
  hJet1_jesUp_ = dirJet.make<TH1D>("h1_jesUp","jesUp",100,0,3);
  hJet1_jesDn_ = dirJet.make<TH1D>("h1_jesDn","jesDn",100,0,3);
  hJet1_jer_   = dirJet.make<TH1D>("h1_jer"  ,"jer"  ,100,0,3);
  hJet1_jerUp_ = dirJet.make<TH1D>("h1_jerUp","jerUp",100,0,3);
  hJet1_jerDn_ = dirJet.make<TH1D>("h1_jerDn","jerDn",100,0,3);
  hJet1_CSV_ = dirJet.make<TH1D>("h1_CSV","CSV",100,-1,1);
  hJet2_pt_    = dirJet.make<TH1D>("h2_pt","pt",500,0,500);
  hJet2_eta_   = dirJet.make<TH1D>("h2_eta","eta",400,-3,3);
  hJet2_phi_   = dirJet.make<TH1D>("h2_phi","phi",100,-4,4);
  hJet2_mass_  = dirJet.make<TH1D>("h2_mass","mass",1000,0,100);
  hJet2_jesUp_ = dirJet.make<TH1D>("h2_jesUp","jesUp",100,0,3);
  hJet2_jesDn_ = dirJet.make<TH1D>("h2_jesDn","jesDn",100,0,3);
  hJet2_jer_   = dirJet.make<TH1D>("h2_jer"  ,"jer"  ,100,0,3);
  hJet2_jerUp_ = dirJet.make<TH1D>("h2_jerUp","jerUp",100,0,3);
  hJet2_jerDn_ = dirJet.make<TH1D>("h2_jerDn","jerDn",100,0,3);
  hJet2_CSV_ = dirJet.make<TH1D>("h2_CSV","CSV",100,-1,1);
  hJet3_pt_    = dirJet.make<TH1D>("h3_pt","pt",500,0,500);
  hJet3_eta_   = dirJet.make<TH1D>("h3_eta","eta",400,-3,3);
  hJet3_phi_   = dirJet.make<TH1D>("h3_phi","phi",100,-4,4);
  hJet3_mass_  = dirJet.make<TH1D>("h3_mass","mass",1000,0,100);
  hJet3_jesUp_ = dirJet.make<TH1D>("h3_jesUp","jesUp",100,0,3);
  hJet3_jesDn_ = dirJet.make<TH1D>("h3_jesDn","jesDn",100,0,3);
  hJet3_jer_   = dirJet.make<TH1D>("h3_jer"  ,"jer"  ,100,0,3);
  hJet3_jerUp_ = dirJet.make<TH1D>("h3_jerUp","jerUp",100,0,3);
  hJet3_jerDn_ = dirJet.make<TH1D>("h3_jerDn","jerDn",100,0,3);
  hJet3_CSV_ = dirJet.make<TH1D>("h3_CSV","CSV",100,-1,1);
  hJet4_pt_    = dirJet.make<TH1D>("h4_pt","pt",500,0,500);
  hJet4_eta_   = dirJet.make<TH1D>("h4_eta","eta",400,-3,3);
  hJet4_phi_   = dirJet.make<TH1D>("h4_phi","phi",100,-4,4);
  hJet4_mass_  = dirJet.make<TH1D>("h4_mass","mass",1000,0,100);
  hJet4_jesUp_ = dirJet.make<TH1D>("h4_jesUp","jesUp",100,0,3);
  hJet4_jesDn_ = dirJet.make<TH1D>("h4_jesDn","jesDn",100,0,3);
  hJet4_jer_   = dirJet.make<TH1D>("h4_jer"  ,"jer"  ,100,0,3);
  hJet4_jerUp_ = dirJet.make<TH1D>("h4_jerUp","jerUp",100,0,3);
  hJet4_jerDn_ = dirJet.make<TH1D>("h4_jerDn","jerDn",100,0,3);
  hJet4_CSV_ = dirJet.make<TH1D>("h4_CSV","CSV",100,-1,1);

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

  const float genWeight = !isMC_ ? 1 : [&](){
    edm::Handle<float> h;
    iEvent.getByToken(genWeightToken_, h);
    return *h;
  }();
  const float puWeight = !isMC_ ? 1 : [&](){
    edm::Handle<float> h;
    iEvent.getByToken(puWeightToken_, h);
    return *h;
  }();
  const double weight = genWeight*puWeight;

  const int nPV = [&](){
    edm::Handle<int> h;
    iEvent.getByToken(nPVToken_, h);
    return *h;
  }();
  hNPV_->Fill(nPV, weight);
  hNPVRaw_->Fill(nPV, genWeight);

  edm::Handle<cat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken_,electronHandle);
  const int nElectron = electronHandle->size();
  hElectron_n_->Fill(nElectron, weight);
  for ( int i=0; i<nElectron; ++i ) {
    const auto& electron = electronHandle->at(i);

    hElectron_pt_   ->Fill( electron.pt() , weight);
    hElectron_eta_  ->Fill( electron.eta() , weight);
    hElectron_phi_  ->Fill( electron.phi() , weight);
    hElectron_mass_ ->Fill( electron.mass() , weight);
    hElectron_chIso_->Fill( electron.chargedHadronIso() , weight);
    hElectron_nhIso_->Fill( electron.neutralHadronIso() , weight);
    hElectron_phIso_->Fill( electron.photonIso() , weight);
    hElectron_puIso_->Fill( electron.puChargedHadronIso() , weight);

    if ( electron.isVeto()      ) hElectron_ids_->Fill(1, weight);
    if ( electron.isMediumMVA() ) hElectron_ids_->Fill(2, weight);
    if ( electron.isTightMVA()  ) hElectron_ids_->Fill(3, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose")  > 0.5 )  hElectron_ids_->Fill(4, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0.5 )  hElectron_ids_->Fill(5, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight")  > 0.5 )  hElectron_ids_->Fill(6, weight);
    if ( electron.isTrigMVAValid() and electron.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.5 )  hElectron_ids_->Fill(7, weight);
    if ( electron.isTrigMVAValid() and electron.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") > 0.5 )  hElectron_ids_->Fill(8, weight);

    hElectron_enDn_->Fill(electron.shiftedEnDown(), weight);
    hElectron_enUp_->Fill(electron.shiftedEnUp(), weight);
  }
  if ( nElectron >= 1 ) {
    const auto& electron = electronHandle->at(0);
    hElectron1_pt_   ->Fill( electron.pt() , weight);
    hElectron1_eta_  ->Fill( electron.eta() , weight);
    hElectron1_phi_  ->Fill( electron.phi() , weight);
    hElectron1_mass_ ->Fill( electron.mass() , weight);
    hElectron1_chIso_->Fill( electron.chargedHadronIso() , weight);
    hElectron1_nhIso_->Fill( electron.neutralHadronIso() , weight);
    hElectron1_phIso_->Fill( electron.photonIso() , weight);
    hElectron1_puIso_->Fill( electron.puChargedHadronIso() , weight);

    if ( electron.isVeto()      ) hElectron1_ids_->Fill(1, weight);
    if ( electron.isMediumMVA() ) hElectron1_ids_->Fill(2, weight);
    if ( electron.isTightMVA()  ) hElectron1_ids_->Fill(3, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose")  > 0.5 )  hElectron1_ids_->Fill(4, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0.5 )  hElectron1_ids_->Fill(5, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight")  > 0.5 )  hElectron1_ids_->Fill(6, weight);
    if ( electron.isTrigMVAValid() and electron.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.5 )  hElectron1_ids_->Fill(7, weight);
    if ( electron.isTrigMVAValid() and electron.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") > 0.5 )  hElectron1_ids_->Fill(8, weight);

    hElectron1_enDn_->Fill(electron.shiftedEnDown(), weight);
    hElectron1_enUp_->Fill(electron.shiftedEnUp(), weight);
  }
  if ( nElectron >= 2 ) {
    const auto& electron = electronHandle->at(1);
    hElectron2_pt_   ->Fill( electron.pt() , weight);
    hElectron2_eta_  ->Fill( electron.eta() , weight);
    hElectron2_phi_  ->Fill( electron.phi() , weight);
    hElectron2_mass_ ->Fill( electron.mass() , weight);
    hElectron2_chIso_->Fill( electron.chargedHadronIso() , weight);
    hElectron2_nhIso_->Fill( electron.neutralHadronIso() , weight);
    hElectron2_phIso_->Fill( electron.photonIso() , weight);
    hElectron2_puIso_->Fill( electron.puChargedHadronIso() , weight);

    if ( electron.isVeto()      ) hElectron2_ids_->Fill(1, weight);
    if ( electron.isMediumMVA() ) hElectron2_ids_->Fill(2, weight);
    if ( electron.isTightMVA()  ) hElectron2_ids_->Fill(3, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose")  > 0.5 )  hElectron2_ids_->Fill(4, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0.5 )  hElectron2_ids_->Fill(5, weight);
    if ( electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight")  > 0.5 )  hElectron2_ids_->Fill(6, weight);
    if ( electron.isTrigMVAValid() and electron.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.5 )  hElectron2_ids_->Fill(7, weight);
    if ( electron.isTrigMVAValid() and electron.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") > 0.5 )  hElectron2_ids_->Fill(8, weight);

    hElectron2_enDn_->Fill(electron.shiftedEnDown(), weight);
    hElectron2_enUp_->Fill(electron.shiftedEnUp(), weight);
  }

  edm::Handle<cat::MuonCollection> muonHandle;
  iEvent.getByToken(muonToken_,muonHandle);
  const int nMuon = muonHandle->size();
  hMuon_n_->Fill(nMuon, weight);
  for ( int i=0; i<nMuon; ++i ) {
    const auto& muon = muonHandle->at(i);

    hMuon_pt_   ->Fill( muon.pt() , weight);
    hMuon_eta_  ->Fill( muon.eta() , weight);
    hMuon_phi_  ->Fill( muon.phi() , weight);
    hMuon_mass_ ->Fill( muon.mass() , weight);
    hMuon_chIso_->Fill( muon.chargedHadronIso() , weight);
    hMuon_nhIso_->Fill( muon.neutralHadronIso() , weight);
    hMuon_phIso_->Fill( muon.photonIso() , weight);
    hMuon_puIso_->Fill( muon.puChargedHadronIso() , weight);

    if ( muon.isLooseMuon() ) hMuon_ids_->Fill(1, weight);
    if ( muon.isTightMuon() ) hMuon_ids_->Fill(2, weight);

    hMuon_enDn_->Fill(muon.shiftedEnDown(), weight);
    hMuon_enUp_->Fill(muon.shiftedEnUp(), weight);
  }
  if ( nMuon >= 1 ) {
    const auto& muon = muonHandle->at(0);
    hMuon1_pt_   ->Fill( muon.pt() , weight);
    hMuon1_eta_  ->Fill( muon.eta() , weight);
    hMuon1_phi_  ->Fill( muon.phi() , weight);
    hMuon1_mass_ ->Fill( muon.mass() , weight);
    hMuon1_chIso_->Fill( muon.chargedHadronIso() , weight);
    hMuon1_nhIso_->Fill( muon.neutralHadronIso() , weight);
    hMuon1_phIso_->Fill( muon.photonIso() , weight);
    hMuon1_puIso_->Fill( muon.puChargedHadronIso() , weight);

    if ( muon.isLooseMuon() ) hMuon1_ids_->Fill(1, weight);
    if ( muon.isTightMuon() ) hMuon1_ids_->Fill(2, weight);

    hMuon1_enDn_->Fill(muon.shiftedEnDown(), weight);
    hMuon1_enUp_->Fill(muon.shiftedEnUp(), weight);
  }
  if ( nMuon >= 2 ) {
    const auto& muon = muonHandle->at(1);
    hMuon2_pt_   ->Fill( muon.pt() , weight);
    hMuon2_eta_  ->Fill( muon.eta() , weight);
    hMuon2_phi_  ->Fill( muon.phi() , weight);
    hMuon2_mass_ ->Fill( muon.mass() , weight);
    hMuon2_chIso_->Fill( muon.chargedHadronIso() , weight);
    hMuon2_nhIso_->Fill( muon.neutralHadronIso() , weight);
    hMuon2_phIso_->Fill( muon.photonIso() , weight);
    hMuon2_puIso_->Fill( muon.puChargedHadronIso() , weight);

    if ( muon.isLooseMuon() ) hMuon2_ids_->Fill(1, weight);
    if ( muon.isTightMuon() ) hMuon2_ids_->Fill(2, weight);

    hMuon2_enDn_->Fill(muon.shiftedEnDown(), weight);
    hMuon2_enUp_->Fill(muon.shiftedEnUp(), weight);
  }

  edm::Handle<cat::PhotonCollection> photonHandle;
  iEvent.getByToken(photonToken_,photonHandle);
  const int nPhoton = photonHandle->size();
  hPhoton_n_->Fill(nPhoton, weight);
  for ( int i=0; i<nPhoton; ++i ) {
    const auto& photon = photonHandle->at(i);

    hPhoton_pt_   ->Fill( photon.pt() , weight);
    hPhoton_eta_  ->Fill( photon.eta() , weight);
    hPhoton_phi_  ->Fill( photon.phi() , weight);
    hPhoton_mass_ ->Fill( photon.mass() , weight);
  }
  if ( nPhoton >= 1 ) {
    const auto& photon = photonHandle->at(0);
    hPhoton1_pt_   ->Fill( photon.pt() , weight);
    hPhoton1_eta_  ->Fill( photon.eta() , weight);
    hPhoton1_phi_  ->Fill( photon.phi() , weight);
    hPhoton1_mass_ ->Fill( photon.mass() , weight);
  }
  if ( nPhoton >= 2 ) {
    const auto& photon = photonHandle->at(1);
    hPhoton2_pt_   ->Fill( photon.pt() , weight);
    hPhoton2_eta_  ->Fill( photon.eta() , weight);
    hPhoton2_phi_  ->Fill( photon.phi() , weight);
    hPhoton2_mass_ ->Fill( photon.mass() , weight);
  }

/*
  edm::Handle<cat::TauCollection> tauHandle;
  iEvent.getByToken(tauToken_,tauHandle);
  const int nTau = tauHandle->size();
  hTau_n_->Fill(nTau, weight);
  for ( int i=0; i<nTau; ++i ) {
    const auto& tau = tauHandle->at(i);

    hTau_pt_   ->Fill( tau.pt() , weight);
    hTau_eta_  ->Fill( tau.eta() , weight);
    hTau_phi_  ->Fill( tau.phi() , weight);
    hTau_mass_ ->Fill( tau.mass() , weight);
  }
  if ( nTau >= 1 ) {
    const auto& tau = tauHandle->at(0);
    hTau1_pt_   ->Fill( tau.pt() , weight);
    hTau1_eta_  ->Fill( tau.eta() , weight);
    hTau1_phi_  ->Fill( tau.phi() , weight);
    hTau1_mass_ ->Fill( tau.mass() , weight);
  }
  if ( nTau >= 2 ) {
    const auto& tau = tauHandle->at(1);
    hTau2_pt_   ->Fill( tau.pt() , weight);
    hTau2_eta_  ->Fill( tau.eta() , weight);
    hTau2_phi_  ->Fill( tau.phi() , weight);
    hTau2_mass_ ->Fill( tau.mass() , weight);
  }
*/

  edm::Handle<cat::JetCollection> jetHandle;
  iEvent.getByToken(jetToken_,jetHandle);
  const int nJet = jetHandle->size();
  hJet_n_->Fill(nJet, weight);
  for ( int i=0; i<nJet; ++i ) {
    const auto& jet = jetHandle->at(i);

    hJet_pt_   ->Fill( jet.pt() , weight);
    hJet_eta_  ->Fill( jet.eta() , weight);
    hJet_phi_  ->Fill( jet.phi() , weight);
    hJet_mass_ ->Fill( jet.mass() , weight);

    hJet_jesUp_->Fill(jet.shiftedEnUp(), weight);
    hJet_jesDn_->Fill(jet.shiftedEnDown(), weight);

    hJet_jer_->Fill(jet.smearedRes(0), weight);
    hJet_jerUp_->Fill(jet.smearedRes(1), weight);
    hJet_jerDn_->Fill(jet.smearedRes(-1), weight);

    hJet_CSV_->Fill(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), weight);
  }
  if ( nJet >= 1 ) {
    const auto& jet = jetHandle->at(0);
    hJet1_pt_   ->Fill( jet.pt() , weight);
    hJet1_eta_  ->Fill( jet.eta() , weight);
    hJet1_phi_  ->Fill( jet.phi() , weight);
    hJet1_mass_ ->Fill( jet.mass() , weight);

    hJet1_jesUp_->Fill(jet.shiftedEnUp(), weight);
    hJet1_jesDn_->Fill(jet.shiftedEnDown(), weight);

    hJet1_jer_->Fill(jet.smearedRes(0), weight);
    hJet1_jerUp_->Fill(jet.smearedRes(1), weight);
    hJet1_jerDn_->Fill(jet.smearedRes(-1), weight);

    hJet_CSV_->Fill(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), weight);
  }
  if ( nJet >= 2 ) {
    const auto& jet = jetHandle->at(1);
    hJet2_pt_   ->Fill( jet.pt() , weight);
    hJet2_eta_  ->Fill( jet.eta() , weight);
    hJet2_phi_  ->Fill( jet.phi() , weight);
    hJet2_mass_ ->Fill( jet.mass() , weight);

    hJet2_jesUp_->Fill(jet.shiftedEnUp(), weight);
    hJet2_jesDn_->Fill(jet.shiftedEnDown(), weight);

    hJet2_jer_->Fill(jet.smearedRes(0), weight);
    hJet2_jerUp_->Fill(jet.smearedRes(1), weight);
    hJet2_jerDn_->Fill(jet.smearedRes(-1), weight);

    hJet_CSV_->Fill(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), weight);
  }
  if ( nJet >= 3 ) {
    const auto& jet = jetHandle->at(2);
    hJet3_pt_   ->Fill( jet.pt() , weight);
    hJet3_eta_  ->Fill( jet.eta() , weight);
    hJet3_phi_  ->Fill( jet.phi() , weight);
    hJet3_mass_ ->Fill( jet.mass() , weight);

    hJet3_jesUp_->Fill(jet.shiftedEnUp(), weight);
    hJet3_jesDn_->Fill(jet.shiftedEnDown(), weight);

    hJet3_jer_->Fill(jet.smearedRes(0), weight);
    hJet3_jerUp_->Fill(jet.smearedRes(1), weight);
    hJet3_jerDn_->Fill(jet.smearedRes(-1), weight);

    hJet3_CSV_->Fill(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), weight);
  }
  if ( nJet >= 4 ) {
    const auto& jet = jetHandle->at(3);
    hJet4_pt_   ->Fill( jet.pt() , weight);
    hJet4_eta_  ->Fill( jet.eta() , weight);
    hJet4_phi_  ->Fill( jet.phi() , weight);
    hJet4_mass_ ->Fill( jet.mass() , weight);

    hJet4_jesUp_->Fill(jet.shiftedEnUp(), weight);
    hJet4_jesDn_->Fill(jet.shiftedEnDown(), weight);

    hJet4_jer_->Fill(jet.smearedRes(0), weight);
    hJet4_jerUp_->Fill(jet.smearedRes(1), weight);
    hJet4_jerDn_->Fill(jet.smearedRes(-1), weight);

    hJet4_CSV_->Fill(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), weight);
  }

  if ( isMC_ ) {
    edm::Handle<reco::GenJetCollection> genJetHandle;
    iEvent.getByToken(genJetToken_,genJetHandle);
    const int nGenJet = genJetHandle->size();
    hGenJet_n_->Fill(nGenJet, weight);
    for ( int i=0; i<nGenJet; ++i ) {
      const auto& genJet = genJetHandle->at(i);

      hGenJet_pt_   ->Fill( genJet.pt() , weight);
      hGenJet_eta_  ->Fill( genJet.eta() , weight);
      hGenJet_phi_  ->Fill( genJet.phi() , weight);
      hGenJet_mass_ ->Fill( genJet.mass() , weight);
    }
    if ( nGenJet >= 1 ) {
      const auto& genJet = genJetHandle->at(0);
      hGenJet1_pt_   ->Fill( genJet.pt() , weight);
      hGenJet1_eta_  ->Fill( genJet.eta() , weight);
      hGenJet1_phi_  ->Fill( genJet.phi() , weight);
      hGenJet1_mass_ ->Fill( genJet.mass() , weight);
    }
    if ( nGenJet >= 2 ) {
      const auto& genJet = genJetHandle->at(1);
      hGenJet2_pt_   ->Fill( genJet.pt() , weight);
      hGenJet2_eta_  ->Fill( genJet.eta() , weight);
      hGenJet2_phi_  ->Fill( genJet.phi() , weight);
      hGenJet2_mass_ ->Fill( genJet.mass() , weight);
    }
    if ( nGenJet >= 3 ) {
      const auto& genJet = genJetHandle->at(2);
      hGenJet3_pt_   ->Fill( genJet.pt() , weight);
      hGenJet3_eta_  ->Fill( genJet.eta() , weight);
      hGenJet3_phi_  ->Fill( genJet.phi() , weight);
      hGenJet3_mass_ ->Fill( genJet.mass() , weight);
    }
    if ( nGenJet >= 4 ) {
      const auto& genJet = genJetHandle->at(3);
      hGenJet4_pt_   ->Fill( genJet.pt() , weight);
      hGenJet4_eta_  ->Fill( genJet.eta() , weight);
      hGenJet4_phi_  ->Fill( genJet.phi() , weight);
      hGenJet4_mass_ ->Fill( genJet.mass() , weight);
    }
  }

  edm::Handle<cat::METCollection> metHandle;
  iEvent.getByToken(metToken_,metHandle);
  if ( !metHandle->empty() ) {
    const auto& met = metHandle->at(0);
    hMET_pt_   ->Fill( met.pt() , weight);
    hMET_phi_  ->Fill( met.phi() , weight);
  }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATHisAnalysis);

