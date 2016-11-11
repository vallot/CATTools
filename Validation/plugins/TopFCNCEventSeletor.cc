#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
//#include "CATTools/DataFormats/interface/SecVertex.h"
//#include "CATTools/DataFormats/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TH1D.h"
#include "TH2F.h"

#define nCutstep 7

using namespace std;

namespace cat {

struct ControlPlotsTopFCNC
{
  typedef TH1D* H1;
  typedef TH2D* H2;

  H1 hCutstep, hCutstepNoweight;

  H1 h_vertex_n[nCutstep];
  H1 h_met_pt[nCutstep], h_met_phi[nCutstep];

  H1 h_leptons_n[nCutstep], h_vetoLeptons_n[nCutstep];
  H1 h_lepton1_pt[nCutstep], h_lepton1_eta[nCutstep], h_lepton1_phi[nCutstep], h_lepton1_q[nCutstep];

  H1 h_electrons_n[nCutstep], h_vetoElectrons_n[nCutstep];
  H1 h_electron1_pt[nCutstep], h_electron1_eta[nCutstep];
  H1 h_electron2_pt[nCutstep], h_electron2_eta[nCutstep];

  H1 h_muons_n[nCutstep], h_vetoMuons_n[nCutstep];
  H1 h_muon1_pt[nCutstep], h_muon1_eta[nCutstep];
  H1 h_muon2_pt[nCutstep], h_muon2_eta[nCutstep];

  H1 h_jets_n[nCutstep], h_jets_pt[nCutstep], h_jets_eta[nCutstep], h_jets_ht[nCutstep];
  H1 h_jet_m[nCutstep][6]; 
  H1 h_jet_pt[nCutstep][6];
  H1 h_jet_eta[nCutstep][6];
  H1 h_jet_phi[nCutstep][6];
  H1 h_jet_btag[nCutstep][6];

  H1 h_bjets_n[nCutstep];
  H1 h_event_st[nCutstep];

  void book(TFileDirectory&& dir)
  {
    const double maxeta = 3;
    const double pi = 3.141592;

    hCutstep = dir.make<TH1D>("cutstep", "cutstep", nCutstep, 1, nCutstep+1);
    hCutstepNoweight = dir.make<TH1D>("cutstepNoweight", "cutstepNoweight", nCutstep, 1, nCutstep+1);

    const char* stepLabels[nCutstep] = {
      "S1 All event", "S2 Good PV", "S3 Trigger",
      "S4 One signal lepton", "S5 Veto other lep", "S6 Veto same lep",
      "S7 nJet3",
    };
    const char* stepNames[nCutstep] = {"step1", "step2", "step3", "step4", "step5", "step6", "step7"};

    for ( int i=0; i<nCutstep; ++i ) {
      hCutstep->GetXaxis()->SetBinLabel(i+1, stepLabels[i]);
      hCutstepNoweight->GetXaxis()->SetBinLabel(i+1, stepLabels[i]);
    }

    std::vector<TFileDirectory> subdirs;

    subdirs.push_back(dir.mkdir(stepNames[0]));
    auto subdir = subdirs.back();
    h_vertex_n[0] = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h_leptons_n[0] = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);

    for ( int i=1; i<=2; ++i ) {
      subdirs.push_back(dir.mkdir(stepNames[i]));
      subdir = subdirs.back();
      h_vertex_n[i] = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
      h_met_pt[i] = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
      h_met_phi[i] = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
      h_leptons_n[i] = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);
      h_vetoLeptons_n[i] = subdir.make<TH1D>("vetoLeptons_n", "vetoLeptons_n", 10, 0, 10);
      h_lepton1_pt[i]  = subdir.make<TH1D>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
      h_lepton1_eta[i] = subdir.make<TH1D>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
      h_lepton1_phi[i] = subdir.make<TH1D>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
      h_lepton1_q[i]   = subdir.make<TH1D>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

      h_jets_n[i] = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
      h_jets_pt[i]  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
      h_jets_eta[i] = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
      h_jets_ht[i] = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);
      h_bjets_n[i] = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);
    }

    for ( int i=3; i<nCutstep; ++i ) {
      subdirs.push_back(dir.mkdir(stepNames[i]));
      subdir = subdirs.back();
      h_vertex_n[i] = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
      h_met_pt[i] = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
      h_met_phi[i] = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
      h_leptons_n[i] = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);

      h_lepton1_pt[i]  = subdir.make<TH1D>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
      h_lepton1_eta[i] = subdir.make<TH1D>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
      h_lepton1_phi[i] = subdir.make<TH1D>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
      h_lepton1_q[i]   = subdir.make<TH1D>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

      h_jets_n[i] = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
      h_jets_pt[i]  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
      h_jets_eta[i] = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
      h_jets_ht[i] = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

      for ( int nJet=0; nJet<6; ++nJet ) {
        const string prefix = Form("jet%d_", nJet+1);
        h_jet_m[i][nJet]   = subdir.make<TH1D>((prefix+"m").c_str(), (prefix+"m").c_str(), 500, 0, 500);
        h_jet_pt[i][nJet]  = subdir.make<TH1D>((prefix+"pt").c_str(), (prefix+"pt").c_str(), 1000, 0, 1000);
        h_jet_eta[i][nJet] = subdir.make<TH1D>((prefix+"eta").c_str(), (prefix+"eta").c_str(), 100, -maxeta, maxeta);
        h_jet_phi[i][nJet] = subdir.make<TH1D>((prefix+"phi").c_str(), (prefix+"phi").c_str(), 100, -pi, pi);
        h_jet_btag[i][nJet] = subdir.make<TH1D>((prefix+"btag").c_str(), (prefix+"btag").c_str(), 100, 0, 1);
      }

      h_bjets_n[i] = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

      h_event_st[i] = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);
    }

    std::vector<int> leptonDebugSteps = {0, 2, 3};
    for ( int step : leptonDebugSteps ) {
      subdir = subdirs[step];
      h_vetoLeptons_n[step] = subdir.make<TH1D>("vetoLeptons_n", "vetoLeptons_n", 10, 0, 10);
      h_electrons_n[step] = subdir.make<TH1D>("electrons_n", "electrons_n", 10, 0, 10);
      h_vetoElectrons_n[step] = subdir.make<TH1D>("vetoElectrons_n", "vetoElectrons_n", 10, 0, 10);
      h_electron1_pt[step]  = subdir.make<TH1D>("electron1_pt", "electron1_pt", 1000, 0, 1000);
      h_electron1_eta[step] = subdir.make<TH1D>("electron1_eta", "electron1_eta", 100, -maxeta, maxeta);
      h_electron2_pt[step]  = subdir.make<TH1D>("electron2_pt", "electron2_pt", 1000, 0, 1000);
      h_electron2_eta[step] = subdir.make<TH1D>("electron2_eta", "electron2_eta", 100, -maxeta, maxeta);
      h_muons_n[step] = subdir.make<TH1D>("muons_n", "muons_n", 10, 0, 10);
      h_vetoMuons_n[step] = subdir.make<TH1D>("vetoMuons_n", "vetoMuons_n", 10, 0, 10);
      h_muon1_pt[step]  = subdir.make<TH1D>("muon1_pt", "muon1_pt", 1000, 0, 1000);
      h_muon1_eta[step] = subdir.make<TH1D>("muon1_eta", "muon1_eta", 100, -maxeta, maxeta);
      h_muon2_pt[step]  = subdir.make<TH1D>("muon2_pt", "muon2_pt", 1000, 0, 1000);
      h_muon2_eta[step] = subdir.make<TH1D>("muon2_eta", "muon2_eta", 100, -maxeta, maxeta);
    }
  };
};

class TopFCNCEventSelector : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TopFCNCEventSelector(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;
  ~TopFCNCEventSelector();

private:
  typedef std::vector<float> vfloat;
  typedef std::vector<double> vdouble;
  edm::EDGetTokenT<float> pileupWeightToken_, genWeightToken_;
  edm::EDGetTokenT<vfloat> genWeightsToken_;
  int genWeightIndex_;

  edm::EDGetTokenT<cat::MuonCollection> muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<cat::METCollection> metToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  edm::EDGetTokenT<int> recoFilterToken_;
  edm::EDGetTokenT<int> trigElToken_, trigMuToken_;
  edm::EDGetTokenT<int> nVertexToken_;

  std::vector<edm::EDGetTokenT<float> > extWeightTokensF_;
  std::vector<edm::EDGetTokenT<double> > extWeightTokensD_;

private:
  TH1D* h_weight, * h_pileupWeight, * h_genWeight;
  ControlPlotsTopFCNC h_el, h_mu;

private:
  double shiftedMuonPt(const cat::Muon& mu) { return mu.pt()+muonScale_*mu.shiftedEn(); }
  double shiftedElectronPt(const cat::Electron& el) { return el.pt()+electronScale_*el.shiftedEn(); }
  double shiftedLepPt(const reco::Candidate& cand)
  {
    auto muonP = dynamic_cast<const cat::Muon*>(&cand);
    auto electronP = dynamic_cast<const cat::Electron*>(&cand);
    if ( muonP ) return shiftedMuonPt(*muonP);
    else if ( electronP ) return shiftedElectronPt(*electronP);
    return cand.pt();
  }
  double shiftedJetPt(const reco::Candidate& cand)
  {
    const auto jet = dynamic_cast<const cat::Jet&>(cand);
    double pt = jet.pt();
    if      ( jetScale_ == +1 ) pt *= jet.shiftedEnUp();
    else if ( jetScale_ == -1 ) pt *= jet.shiftedEnDown();

    if ( isMC_ and !isSkipJER_ ) pt *= jet.smearedRes(jetResol_);

    return pt;
  }

  bool isGoodMuon(const cat::Muon& mu)
  {
    if ( std::abs(mu.eta()) > 2.1 ) return false;
    if ( std::isnan(mu.pt()) or shiftedMuonPt(mu) < 27 ) return false;

    if ( mu.relIso(0.4) > 0.15 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const cat::Electron& el)
  {
    if ( std::abs(el.eta()) > 2.1 ) return false;
    if ( std::isnan(el.pt()) or shiftedElectronPt(el) < 35 ) return false;

    if ( isMVAElectronSel_ and !el.isTrigMVAValid() ) return false;

    //if ( el.relIso(0.3) >= 0.11 ) return false;
    if ( !el.electronID(elIdName_) ) return false;
    //if ( !el.isPF() or !el.passConversionVeto() ) return false;
    //const double scEta = std::abs(el.scEta());
    //if ( isEcalCrackVeto_ and scEta > 1.4442 and scEta < 1.566 ) return false;
    //const double d0 = std::abs(el.dxy()), dz = std::abs(el.vz());
    //if      ( scEta <= 1.479 and (d0 > 0.05 or dz > 0.1) ) return false;
    //else if ( scEta >  1.479 and (d0 > 0.10 or dz > 0.2) ) return false;
    return true;
  }
  bool isVetoMuon(const cat::Muon& mu)
  {
    if ( std::abs(mu.eta()) > 2.4 ) return false;
    if ( std::isnan(mu.pt()) or shiftedMuonPt(mu) < 10 ) return false;

    if ( !mu.isLooseMuon() ) return false;
    if ( mu.relIso(0.4) > 0.25 ) return false;
    return true;
  }
  bool isVetoElectron(const cat::Electron& el)
  {
    if ( std::abs(el.eta()) > 2.4 ) return false;
    if ( std::isnan(el.pt()) or shiftedElectronPt(el) < 10 ) return false;
    if ( !el.electronID(elVetoIdName_) ) return false;
    //const double scEta = std::abs(el.scEta());
    //const double d0 = std::abs(el.dxy()), dz = std::abs(el.vz());
    //if      ( scEta <= 1.479 and (d0 > 0.05 or dz > 0.1) ) return false;
    //else if ( scEta >  1.479 and (d0 > 0.10 or dz > 0.2) ) return false;
    return true;
  }
  bool isBjet(const cat::Jet& jet)
  {
    const double bTag = jet.bDiscriminator(bTagName_);
    if      ( bTagWP_ == BTagWP::CSVL ) return bTag > WP_BTAG_CSVv2L;
    else if ( bTagWP_ == BTagWP::CSVM ) return bTag > WP_BTAG_CSVv2M;
    else if ( bTagWP_ == BTagWP::CSVT ) return bTag > WP_BTAG_CSVv2T;
    return false;
  }

private:
  typedef reco::Candidate::LorentzVector LV;

  // Energy scales
  int muonScale_, electronScale_, jetScale_, jetResol_;
  bool isSkipJER_; // Do not apply JER, needed to remove randomness during the Synchronization

  // Efficiency SF
  ScaleFactorEvaluator muonSF_, electronSF_;
  double muonSFShift_, electronSFShift_;
  double trigSFShift_;

  bool isMC_;
  bool isIgnoreTrig_; // Accept event even if it does not pass HLT. Needed for synchronization
  const int applyFilterAt_;

  // ID variables
  bool isEcalCrackVeto_, isMVAElectronSel_;
  std::string bTagName_;
  std::string elIdName_, elVetoIdName_;
  enum class BTagWP { CSVL, CSVM, CSVT } bTagWP_;

};

}

using namespace cat;

TopFCNCEventSelector::TopFCNCEventSelector(const edm::ParameterSet& pset):
  isMC_(pset.getParameter<bool>("isMC")),
  applyFilterAt_(pset.getParameter<int>("applyFilterAt"))
{
  const auto muonSet = pset.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  muonScale_ = muonSet.getParameter<int>("scaleDirection");
  if ( isMC_ ) {
    const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("efficiencySF");
    // FIXME : for muons, eta bins are folded - always double check this with cfg
    muonSF_.set(muonSFSet.getParameter<vdouble>("pt_bins"),
                muonSFSet.getParameter<vdouble>("abseta_bins"),
                muonSFSet.getParameter<vdouble>("values"),
                muonSFSet.getParameter<vdouble>("errors"));
    muonSFShift_ = muonSet.getParameter<int>("efficiencySFDirection");
  }

  const auto electronSet = pset.getParameter<edm::ParameterSet>("electron");
  electronToken_ = consumes<cat::ElectronCollection>(electronSet.getParameter<edm::InputTag>("src"));
  elIdName_ = electronSet.getParameter<string>("idName");
  elVetoIdName_ = electronSet.getParameter<string>("vetoIdName");
  electronScale_ = electronSet.getParameter<int>("scaleDirection");
  if ( isMC_ ) {
    const auto electronSFSet = electronSet.getParameter<edm::ParameterSet>("efficiencySF");
    // FIXME : for electrons, eta bins are NOT folded - always double check this with cfg
    electronSF_.set(electronSFSet.getParameter<vdouble>("pt_bins"),
                    electronSFSet.getParameter<vdouble>("abseta_bins"),
                    electronSFSet.getParameter<vdouble>("values"),
                    electronSFSet.getParameter<vdouble>("errors"));
    electronSFShift_ = electronSet.getParameter<int>("efficiencySFDirection");
  }
  isEcalCrackVeto_ = isMVAElectronSel_ = false;
  if ( elIdName_.substr(0,3) == "mva" ) {
    isMVAElectronSel_ = true;
  }
  else {
    isEcalCrackVeto_ = electronSet.getParameter<bool>("applyEcalCrackVeto");
  }

  const auto jetSet = pset.getParameter<edm::ParameterSet>("jet");
  jetToken_ = consumes<cat::JetCollection>(jetSet.getParameter<edm::InputTag>("src"));
  jetScale_ = jetSet.getParameter<int>("scaleDirection");
  jetResol_ = jetSet.getParameter<int>("resolDirection");
  bTagName_ = jetSet.getParameter<string>("bTagName");
  const auto bTagWPStr = jetSet.getParameter<string>("bTagWP");
  if      ( bTagWPStr == "CSVL" ) bTagWP_ = BTagWP::CSVL;
  else if ( bTagWPStr == "CSVM" ) bTagWP_ = BTagWP::CSVM;
  else if ( bTagWPStr == "CSVT" ) bTagWP_ = BTagWP::CSVT;
  else edm::LogError("TopFCNCEventSelector") << "Wrong bTagWP parameter " << bTagWPStr;
  isSkipJER_ = jetSet.getParameter<bool>("skipJER");

  const auto metSet = pset.getParameter<edm::ParameterSet>("met");
  metToken_ = consumes<cat::METCollection>(metSet.getParameter<edm::InputTag>("src"));

  const auto vertexSet = pset.getParameter<edm::ParameterSet>("vertex");
  nVertexToken_ = consumes<int>(vertexSet.getParameter<edm::InputTag>("nVertex"));
  vertexToken_ = consumes<reco::VertexCollection>(vertexSet.getParameter<edm::InputTag>("src"));
  pileupWeightToken_ = consumes<float>(vertexSet.getParameter<edm::InputTag>("pileupWeight"));

  const auto filterSet = pset.getParameter<edm::ParameterSet>("filters");
  recoFilterToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("filterRECO"));
  trigElToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigEL"));
  trigMuToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigMU"));
  isIgnoreTrig_ = filterSet.getParameter<bool>("ignoreTrig");
  trigSFShift_ = filterSet.getParameter<int>("efficiencySFDirection");

  if ( isMC_ ) {
    const auto genWeightSet = pset.getParameter<edm::ParameterSet>("genWeight");

    genWeightIndex_ = genWeightSet.getParameter<int>("index");
    if ( genWeightIndex_ < 0 ) genWeightToken_ = consumes<float>(genWeightSet.getParameter<edm::InputTag>("src"));
    else genWeightsToken_ = consumes<vfloat>(genWeightSet.getParameter<edm::InputTag>("src"));
  }

  // Other weights
  const auto extWeightLabels = pset.getParameter<std::vector<edm::InputTag> >("extWeights");
  for ( auto x : extWeightLabels ) {
    extWeightTokensF_.push_back(consumes<float>(x));
    extWeightTokensD_.push_back(consumes<double>(x));
  }

  // Fill histograms, etc
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  auto doverall = fs->mkdir("overall", "overall");
  h_weight = doverall.make<TH1D>("weight", "weight", 200, -10, 10);
  if ( isMC_ ) {
    h_genWeight = doverall.make<TH1D>("genWeight", "genWeight", 200, -10, 10);
    h_pileupWeight = doverall.make<TH1D>("pileupWeight", "pileupWeight", 200, -10, 10);
  }

  h_el.book(fs->mkdir("el"));
  h_mu.book(fs->mkdir("mu"));

  produces<int>("channel");
  produces<float>("weight");
  produces<float>("met");
  produces<float>("metphi");
  produces<std::vector<cat::Lepton> >("leptons");
  produces<std::vector<cat::Jet> >("jets");
}

bool TopFCNCEventSelector::filter(edm::Event& event, const edm::EventSetup&)
{
  if ( event.isRealData() ) isMC_ = false;

  // Get physics objects
  edm::Handle<cat::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<cat::ElectronCollection> electronHandle;
  event.getByToken(electronToken_, electronHandle);

  edm::Handle<cat::JetCollection> jetHandle;
  event.getByToken(jetToken_, jetHandle);

  edm::Handle<cat::METCollection> metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& metP4 = metHandle->at(0).p4();
  double metDpx = 0, metDpy = 0;

  //edm::Handle<int> nVertexHandle;
  //event.getByToken(nVertexToken_, nVertexHandle);
  //const int nVertex = *nVertexHandle;
  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  const int nVertex = [&](){
    int n = 0;
    for ( auto& v : *vertexHandle ) {
      if ( v.isFake() or v.ndof() <= 4 ) continue;
      //if ( std::abs(v.position().rho()) >= 2 or std::abs(v.z()) >= 24 ) continue;     
      if ( std::abs(v.position().rho()) >= 2 or std::abs(v.z()) >= 20 ) continue;     
      ++n;
    }
    return n;
  }();

  std::auto_ptr<std::vector<cat::Lepton> > out_leptons(new std::vector<cat::Lepton>());
  std::auto_ptr<std::vector<cat::Jet> > out_jets(new std::vector<cat::Jet>());

  // Compute event weight - from generator, pileup, etc
  double weight = 1.0;
  if ( isMC_ ) {
    float genWeight = 1.;
    edm::Handle<float> fHandle;
    edm::Handle<vfloat> vfHandle;

    if ( genWeightIndex_ < 0 ) {
      event.getByToken(genWeightToken_, fHandle);
      genWeight = *fHandle;
    }
    else {
      event.getByToken(genWeightsToken_, vfHandle);
      genWeight = vfHandle->at(genWeightIndex_);
    }

    event.getByToken(pileupWeightToken_, fHandle);
    const float pileupWeight = *fHandle;

    h_genWeight->Fill(genWeight);
    h_pileupWeight->Fill(pileupWeight);
    weight *= genWeight*pileupWeight;
    // NOTE: weight value to be multiplied by lepton SF, etc.
  }

  // Apply all other weights
  for ( auto t : extWeightTokensF_ ) {
    edm::Handle<float> h;
    if ( event.getByToken(t, h) ) weight *= *h;
  }
  for ( auto t : extWeightTokensD_ ) {
    edm::Handle<double> h;
    if ( event.getByToken(t, h) ) weight *= *h;
  }

  // Get event filters and triggers
  edm::Handle<int> trigHandle;
  event.getByToken(recoFilterToken_, trigHandle);
  //const int isRECOFilterOK = *trigHandle;

  event.getByToken(trigElToken_, trigHandle);
  const int isTrigEl = *trigHandle;
  event.getByToken(trigMuToken_, trigHandle);
  const int isTrigMu = *trigHandle;

  // Select good leptons
  double leptons_st = 0;
  cat::MuonCollection selMuons, vetoMuons;
  for ( int i=0, n=muonHandle->size(); i<n; ++i ) {
    auto& p = muonHandle->at(i);
    const double pt = shiftedMuonPt(p);
    const double scale = pt/p.pt();

    cat::Muon lep(p);
    lep.setP4(p.p4()*scale);
    const bool isGood = isGoodMuon(p); // note: pt scale is done in the function
    const bool isVeto = isVetoMuon(p); // note: pt scale is done in the function
    if ( isGood ) selMuons.push_back(lep);
    else if ( isVeto ) vetoMuons.push_back(lep);
    else continue;

    leptons_st += pt;
    metDpx += lep.px()-p.px();
    metDpy += lep.py()-p.py();
  }
  cat::ElectronCollection selElectrons, vetoElectrons;
  for ( int i=0, n=electronHandle->size(); i<n; ++i ) {
    auto& p = electronHandle->at(i);
    const double pt = shiftedElectronPt(p);
    const double scale = pt/p.pt();

    cat::Electron lep(p);
    lep.setP4(p.p4()*scale);
    const bool isGood = isGoodElectron(p); // note: pt scale is done in the function
    const bool isVeto = isVetoElectron(p); // note: pt scale is done in the function
    if ( isGood ) selElectrons.push_back(lep);
    else if ( isVeto ) vetoElectrons.push_back(lep);
    else continue;

    leptons_st += pt;
    metDpx += lep.px()-p.px();
    metDpy += lep.py()-p.py();
  }
  std::vector<const cat::Lepton*> selLeptons;
  for ( auto& x : selMuons ) selLeptons.push_back(&x);
  for ( auto& x : selElectrons ) selLeptons.push_back(&x);
  std::sort(selLeptons.begin(), selLeptons.end(),
            [&](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
  // Copy selLeptons to out_leptons
  if ( !selLeptons.empty() ) out_leptons->push_back(*selLeptons.at(0));
  const int leptons_n = selLeptons.size();
  const cat::Lepton* lepton1 = 0;
  int channel = 0;
  if ( leptons_n >= 1 ) {
    // Set lepton1
    lepton1 = selLeptons.at(0);

    const int pdgId1 = std::abs(lepton1->pdgId());
    // Determine channel
    channel = abs(pdgId1);

    // Apply lepton SF
    if ( channel == 11 ) {
      const auto e1 = dynamic_cast<const cat::Electron*>(lepton1);
      const double w1 = electronSF_(lepton1->pt(), std::abs(e1->scEta()), electronSFShift_);
      weight *= w1;
      if ( !isIgnoreTrig_ ) weight *= isTrigEl;// * computeTrigSF(*lepton1, trigSFShift_);
    }
    else if ( channel == 13 ) {
      const double w1 = muonSF_(lepton1->pt(), std::abs(lepton1->eta()), muonSFShift_);
      weight *= w1;
      if ( !isIgnoreTrig_ ) weight *= isTrigMu;// * computeTrigSF(*lepton1, trigSFShift_);
    }
    else edm::LogError("TopFCNCEventSelector") << "Strange event with nLepton >=2 but not falling info ee,mumu,emu category";
  }

  // Select good jets
  int bjets_n = 0;
  double jets_ht = 0;
  for ( int i=0, n=jetHandle->size(); i<n; ++i ) {
    auto& p = jetHandle->at(i);
    if ( std::abs(p.eta()) > 2.4 ) continue;
    if ( !p.LooseId() ) continue;

    const double pt = shiftedJetPt(p);
    const double scale = pt/p.pt();
    cat::Jet jet(p);
    jet.setP4(scale*p.p4());

    metDpx += jet.px()-p.px();
    metDpy += jet.py()-p.py();
    if ( pt < 30 ) continue;

    if ( leptons_n >= 1 and deltaR(jet.p4(), out_leptons->at(0).p4()) < 0.4 ) continue;

    out_jets->push_back(jet);
    jets_ht += pt;
    if ( isBjet(p) ) ++bjets_n;
  }
  const int jets_n = out_jets->size();
  std::sort(out_jets->begin(), out_jets->end(),
            [&](const cat::Jet& a, const cat::Jet& b){return a.pt() > b.pt();});

  // Update & calculate met
  const double met_pt = hypot(metP4.px()-metDpx, metP4.py()-metDpy);
  const double met_phi = atan2(metP4.px()-metDpx, metP4.py()-metDpy);

  // Check cut steps and fill some histograms
  int cutstep = 1;
  switch(1) default: { // C++ trick to avoid nested if-statements.
    // Fill all events
    h_weight->Fill(weight);

    h_el.hCutstep->Fill(cutstep, weight);
    h_el.hCutstepNoweight->Fill(cutstep);

    h_el.h_vertex_n[0]->Fill(nVertex, weight);
    h_el.h_leptons_n[0]->Fill(leptons_n, weight);
    h_el.h_vetoLeptons_n[0]->Fill(vetoMuons.size()+vetoElectrons.size(), weight);
    h_el.h_electrons_n[0]->Fill(selElectrons.size(), weight);
    h_el.h_vetoElectrons_n[0]->Fill(vetoElectrons.size(), weight);
    h_el.h_muons_n[0]->Fill(selMuons.size(), weight);
    h_el.h_vetoMuons_n[0]->Fill(vetoMuons.size(), weight);

    if ( selElectrons.size() > 0 ) {
      const auto p = selElectrons.at(0);
      h_el.h_electron1_pt[0]->Fill(p.pt(), weight);
      h_el.h_electron1_eta[0]->Fill(p.eta(), weight);
    }
    if ( selElectrons.size() > 1 ) {
      const auto p = selElectrons.at(1);
      h_el.h_electron2_pt[0]->Fill(p.pt(), weight);
      h_el.h_electron2_eta[0]->Fill(p.eta(), weight);
    }
    if ( selMuons.size() > 0 ) {
      const auto p = selMuons.at(0);
      h_el.h_muon1_pt[0]->Fill(p.pt(), weight);
      h_el.h_muon1_eta[0]->Fill(p.eta(), weight);
    }
    if ( selMuons.size() > 1 ) {
      const auto p = selMuons.at(1);
      h_el.h_muon2_pt[0]->Fill(p.pt(), weight);
      h_el.h_muon2_eta[0]->Fill(p.eta(), weight);
    }

    h_mu.hCutstep->Fill(cutstep, weight);
    h_mu.hCutstepNoweight->Fill(cutstep);
    h_mu.h_vertex_n[0]->Fill(nVertex, weight);

    h_mu.h_vetoLeptons_n[0]->Fill(vetoMuons.size()+vetoElectrons.size(), weight);
    h_mu.h_electrons_n[0]->Fill(selElectrons.size(), weight);
    h_mu.h_vetoElectrons_n[0]->Fill(vetoElectrons.size(), weight);
    h_mu.h_muons_n[0]->Fill(selMuons.size(), weight);
    h_mu.h_vetoMuons_n[0]->Fill(vetoMuons.size(), weight);

    if ( selElectrons.size() > 0 ) {
      const auto p = selElectrons.at(0);
      h_mu.h_electron1_pt[0]->Fill(p.pt(), weight);
      h_mu.h_electron1_eta[0]->Fill(p.eta(), weight);
    }
    if ( selElectrons.size() > 1 ) {
      const auto p = selElectrons.at(1);
      h_mu.h_electron2_pt[0]->Fill(p.pt(), weight);
      h_mu.h_electron2_eta[0]->Fill(p.eta(), weight);
    }
    if ( selMuons.size() > 0 ) {
      const auto p = selMuons.at(0);
      h_mu.h_muon1_pt[0]->Fill(p.pt(), weight);
      h_mu.h_muon1_eta[0]->Fill(p.eta(), weight);
    }
    if ( selMuons.size() > 1 ) {
      const auto p = selMuons.at(1);
      h_mu.h_muon2_pt[0]->Fill(p.pt(), weight);
      h_mu.h_muon2_eta[0]->Fill(p.eta(), weight);
    }

    if ( nVertex <= 0 ) break;
    cutstep = 2;

    h_el.hCutstep->Fill(cutstep, weight);
    h_el.hCutstepNoweight->Fill(cutstep);
    h_el.h_vertex_n[1]->Fill(nVertex, weight);
    h_el.h_met_pt[1]->Fill(met_pt, weight);
    h_el.h_met_phi[1]->Fill(met_phi, weight);
    h_el.h_leptons_n[1]->Fill(leptons_n, weight);

    if ( leptons_n >= 1 ) {
      const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();
      h_el.h_lepton1_pt[1]->Fill(lepton1P4.pt(), weight);
      h_el.h_lepton1_eta[1]->Fill(lepton1->eta(), weight);
      h_el.h_lepton1_phi[1]->Fill(lepton1->phi(), weight);
      h_el.h_lepton1_q[1]->Fill(lepton1->charge(), weight);
    }

    h_el.h_jets_n[1]->Fill(jets_n, weight);
    h_el.h_bjets_n[1]->Fill(bjets_n, weight);
    h_el.h_jets_ht[1]->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_el.h_jets_pt[1]->Fill(jet.pt(), weight);
      h_el.h_jets_eta[1]->Fill(jet.eta(), weight);
    }

    h_mu.hCutstep->Fill(cutstep, weight);
    h_mu.hCutstepNoweight->Fill(cutstep);
    h_mu.h_vertex_n[1]->Fill(nVertex, weight);
    h_mu.h_met_pt[1]->Fill(met_pt, weight);
    h_mu.h_met_phi[1]->Fill(met_phi, weight);
    h_mu.h_leptons_n[1]->Fill(leptons_n, weight);
    if ( leptons_n >= 1 ) {
      const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();
      h_mu.h_lepton1_pt[1]->Fill(lepton1P4.pt(), weight);
      h_mu.h_lepton1_eta[1]->Fill(lepton1->eta(), weight);
      h_mu.h_lepton1_phi[1]->Fill(lepton1->phi(), weight);
      h_mu.h_lepton1_q[1]->Fill(lepton1->charge(), weight);
    }
    h_mu.h_jets_n[1]->Fill(jets_n, weight);
    h_mu.h_bjets_n[1]->Fill(bjets_n, weight);
    h_mu.h_jets_ht[1]->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_mu.h_jets_pt[1]->Fill(jet.pt(), weight);
      h_mu.h_jets_eta[1]->Fill(jet.eta(), weight);
    }

    // El channel Cutstep 0b with trigger requirements
    int cutstep_el = cutstep, cutstep_mu = cutstep;
    if ( isIgnoreTrig_ or isTrigEl ) {
      cutstep_el = 3;
      h_el.hCutstep->Fill(cutstep_el, weight);
      h_el.hCutstepNoweight->Fill(cutstep_el);
      h_el.h_vertex_n[2]->Fill(nVertex, weight);
      h_el.h_met_pt[2]->Fill(met_pt, weight);
      h_el.h_met_phi[2]->Fill(met_phi, weight);
      h_el.h_leptons_n[2]->Fill(leptons_n, weight);

      if ( leptons_n >= 1 ) {
        const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();
        h_el.h_lepton1_pt[2]->Fill(lepton1P4.pt(), weight);
        h_el.h_lepton1_eta[2]->Fill(lepton1->eta(), weight);
        h_el.h_lepton1_phi[2]->Fill(lepton1->phi(), weight);
        h_el.h_lepton1_q[2]->Fill(lepton1->charge(), weight);
      }

      h_el.h_vetoLeptons_n[2]->Fill(vetoMuons.size()+vetoElectrons.size(), weight);
      h_el.h_electrons_n[2]->Fill(selElectrons.size(), weight);
      h_el.h_vetoElectrons_n[2]->Fill(vetoElectrons.size(), weight);
      h_el.h_muons_n[2]->Fill(selMuons.size(), weight);
      h_el.h_vetoMuons_n[2]->Fill(vetoMuons.size(), weight);

      if ( selElectrons.size() > 0 ) {
        const auto p = selElectrons.at(0);
        h_el.h_electron1_pt[2]->Fill(p.pt(), weight);
        h_el.h_electron1_eta[2]->Fill(p.eta(), weight);
      }
      if ( selElectrons.size() > 1 ) {
        const auto p = selElectrons.at(1);
        h_el.h_electron2_pt[2]->Fill(p.pt(), weight);
        h_el.h_electron2_eta[2]->Fill(p.eta(), weight);
      }
      if ( selMuons.size() > 0 ) {
        const auto p = selMuons.at(0);
        h_el.h_muon1_pt[2]->Fill(p.pt(), weight);
        h_el.h_muon1_eta[2]->Fill(p.eta(), weight);
      }
      if ( selMuons.size() > 1 ) {
        const auto p = selMuons.at(1);
        h_el.h_muon2_pt[2]->Fill(p.pt(), weight);
        h_el.h_muon2_eta[2]->Fill(p.eta(), weight);
      }

      h_el.h_jets_n[2]->Fill(jets_n, weight);
      h_el.h_bjets_n[2]->Fill(bjets_n, weight);
      h_el.h_jets_ht[2]->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h_jets_pt[2]->Fill(jet.pt(), weight);
        h_el.h_jets_eta[2]->Fill(jet.eta(), weight);
      }

    }
    // Mu channel Cutstep 0b with trigger requirements
    if ( isIgnoreTrig_ or isTrigMu ) {
      cutstep_mu = 3;
      h_mu.hCutstep->Fill(cutstep_mu, weight);
      h_mu.hCutstepNoweight->Fill(cutstep_mu);
      h_mu.h_vertex_n[2]->Fill(nVertex, weight);
      h_mu.h_met_pt[2]->Fill(met_pt, weight);
      h_mu.h_met_phi[2]->Fill(met_phi, weight);
      h_mu.h_leptons_n[2]->Fill(leptons_n, weight);

      h_mu.h_vetoLeptons_n[2]->Fill(vetoMuons.size()+vetoElectrons.size(), weight);
      h_mu.h_electrons_n[2]->Fill(selElectrons.size(), weight);
      h_mu.h_vetoElectrons_n[2]->Fill(vetoElectrons.size(), weight);
      h_mu.h_muons_n[2]->Fill(selMuons.size(), weight);
      h_mu.h_vetoMuons_n[2]->Fill(vetoMuons.size(), weight);

      if ( selElectrons.size() > 0 ) {
        const auto p = selElectrons.at(0);
        h_mu.h_electron1_pt[2]->Fill(p.pt(), weight);
        h_mu.h_electron1_eta[2]->Fill(p.eta(), weight);
      }
      if ( selElectrons.size() > 1 ) {
        const auto p = selElectrons.at(1);
        h_mu.h_electron2_pt[2]->Fill(p.pt(), weight);
        h_mu.h_electron2_eta[2]->Fill(p.eta(), weight);
      }
      if ( selMuons.size() > 0 ) {
        const auto p = selMuons.at(0);
        h_mu.h_muon1_pt[2]->Fill(p.pt(), weight);
        h_mu.h_muon1_eta[2]->Fill(p.eta(), weight);
      }
      if ( selMuons.size() > 1 ) {
        const auto p = selMuons.at(1);
        h_mu.h_muon2_pt[2]->Fill(p.pt(), weight);
        h_mu.h_muon2_eta[2]->Fill(p.eta(), weight);
      }

      if ( leptons_n >= 1 ) {
        const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();
        h_mu.h_lepton1_pt[2]->Fill(lepton1P4.pt(), weight);
        h_mu.h_lepton1_eta[2]->Fill(lepton1->eta(), weight);
        h_mu.h_lepton1_phi[2]->Fill(lepton1->phi(), weight);
        h_mu.h_lepton1_q[2]->Fill(lepton1->charge(), weight);
      }

      h_mu.h_jets_n[2]->Fill(jets_n, weight);
      h_mu.h_bjets_n[2]->Fill(bjets_n, weight);
      h_mu.h_jets_ht[2]->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h_jets_pt[2]->Fill(jet.pt(), weight);
        h_mu.h_jets_eta[2]->Fill(jet.eta(), weight);
      }
    }

    // Apply up to cut step3
    if      ( channel == 11 and cutstep_el == 3 ) cutstep = 3;
    else if ( channel == 13 and cutstep_mu == 3 ) cutstep = 3;
    else break;

    // Now the next steps are rather clear, start from lepton selection

    // Step4 exactly one signal lepton
    if ( leptons_n != 1 ) break;
    cutstep = 4;

    // Step5 veto any additional lepton in different flavour
    if ( (channel == 11 and !vetoMuons.empty()) or
        (channel == 13 and !vetoElectrons.empty()) ) break;
    cutstep = 5;

    // Step6 veto any additional lepton in same flavour
    if ( (channel == 11 and !vetoElectrons.empty()) or
        (channel == 13 and !vetoMuons.empty()) ) break;
    cutstep = 6;

    // Step7 Minimal jet multiplicity
    if ( jets_n < 3 ) break;
    cutstep = 7;
  }

  // Cut step is ready. Now proceed to fill histograms from step 4 one lepton
  if ( cutstep >= 4 ) {
    auto& h = channel == 11 ? h_el : h_mu;
    // Start from the step1 (is [3] in the array)
    // lepton1 should exist from step1
    const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();

    for ( int icutstep=4; icutstep<=cutstep; ++icutstep ) {
      const int i = icutstep-1; // Fill phys object histograms from one lepton step, array[3]
      h.hCutstep->Fill(icutstep, weight);
      h.hCutstepNoweight->Fill(icutstep);
      h.h_vertex_n[i]->Fill(nVertex, weight);
      h.h_met_pt[i]->Fill(met_pt, weight);
      h.h_met_phi[i]->Fill(met_phi, weight);
      h.h_leptons_n[i]->Fill(leptons_n, weight);
      h.h_lepton1_pt[i]->Fill(lepton1P4.pt(), weight);
      h.h_lepton1_eta[i]->Fill(lepton1->eta(), weight);
      h.h_lepton1_phi[i]->Fill(lepton1->phi(), weight);
      h.h_lepton1_q[i]->Fill(lepton1->charge(), weight);
      h.h_jets_n[i]->Fill(jets_n, weight);
      h.h_jets_ht[i]->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h.h_jets_pt[i]->Fill(jet.pt(), weight);
        h.h_jets_eta[i]->Fill(jet.eta(), weight);
      }
      for ( int j=0, n=std::min(6, jets_n); j<n; ++j ) {
        const auto& jet = out_jets->at(j);
        h.h_jet_m[i][j]->Fill(jet.mass(), weight);
        h.h_jet_pt[i][j]->Fill(jet.pt(), weight);
        h.h_jet_eta[i][j]->Fill(jet.eta(), weight);
        h.h_jet_phi[i][j]->Fill(jet.phi(), weight);
        h.h_jet_btag[i][j]->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h.h_bjets_n[i]->Fill(bjets_n, weight);
      h.h_event_st[i]->Fill(leptons_st+jets_ht+met_pt, weight);

      if ( i == 0 or i == 2 or i == 3 ) {
        h.h_vetoLeptons_n[i]->Fill(vetoMuons.size()+vetoElectrons.size(), weight);
        h.h_electrons_n[i]->Fill(selElectrons.size(), weight);
        h.h_vetoElectrons_n[i]->Fill(vetoElectrons.size(), weight);
        h.h_vetoElectrons_n[i]->Fill(selMuons.size(), weight);
        h.h_vetoMuons_n[i]->Fill(vetoMuons.size(), weight);

        if ( selElectrons.size() > 0 ) {
          const auto p = selElectrons.at(0);
          h.h_electron1_pt[i]->Fill(p.pt(), weight);
          h.h_electron1_eta[i]->Fill(p.eta(), weight);
        }
        if ( selElectrons.size() > 1 ) {
          const auto p = selElectrons.at(1);
          h.h_electron2_pt[i]->Fill(p.pt(), weight);
          h.h_electron2_eta[i]->Fill(p.eta(), weight);
        }
        if ( selMuons.size() > 0 ) {
          const auto p = selMuons.at(0);
          h.h_muon1_pt[i]->Fill(p.pt(), weight);
          h.h_muon1_eta[i]->Fill(p.eta(), weight);
        }
        if ( selMuons.size() > 1 ) {
          const auto p = selMuons.at(1);
          h.h_muon2_pt[i]->Fill(p.pt(), weight);
          h.h_muon2_eta[i]->Fill(p.eta(), weight);
        }

        if ( leptons_n >= 1 ) {
          const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();
          h.h_lepton1_pt[i]->Fill(lepton1P4.pt(), weight);
          h.h_lepton1_eta[i]->Fill(lepton1->eta(), weight);
          h.h_lepton1_phi[i]->Fill(lepton1->phi(), weight);
          h.h_lepton1_q[i]->Fill(lepton1->charge(), weight);
        }
      }

    }
  }

  event.put(std::auto_ptr<int>(new int((int)channel)), "channel");
  event.put(std::auto_ptr<float>(new float(weight)), "weight");
  event.put(std::auto_ptr<float>(new float(metP4.pt())), "met");
  event.put(std::auto_ptr<float>(new float(metP4.phi())), "metphi");
  event.put(out_leptons, "leptons");
  event.put(out_jets, "jets");

  // Apply filter at the given step.
  if ( cutstep >= applyFilterAt_ ) return true;

  return false;
}

TopFCNCEventSelector::~TopFCNCEventSelector()
{
  if ( h_el.hCutstepNoweight ) {
    cout << "---- cut flows without weight ----\n";
    cout << "Step\tel\tmu\n";
    const int n = h_el.hCutstepNoweight->GetNbinsX();
    for ( int i=1; i<=n; ++i ) {
      const string name(h_el.hCutstepNoweight->GetXaxis()->GetBinLabel(i));
      if ( name.empty() ) break;
      cout << name.substr(0, name.find(' '));
      cout << '\t' << h_el.hCutstepNoweight->GetBinContent(i);
      cout << '\t' << h_mu.hCutstepNoweight->GetBinContent(i) << '\n';
    }
    cout << "-----------------------------------\n";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopFCNCEventSelector);

