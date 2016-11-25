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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"

#include <fstream>

#define nCutsteps 7

using namespace std;

class TopFCNCEventSelector : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TopFCNCEventSelector(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;
  ~TopFCNCEventSelector();

private:
  typedef std::vector<float> vfloat;
  typedef std::vector<double> vdouble;
  int channel_;
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
  double shiftedMuonScale(const cat::Muon& mu) { return 1+muonScale_*mu.shiftedEn()/mu.pt(); }
  double shiftedElectronScale(const cat::Electron& el) { return 1+electronScale_*el.shiftedEn()/el.pt(); }
  double shiftedJetScale(const reco::Candidate& cand)
  {
    const auto jet = dynamic_cast<const cat::Jet&>(cand);
    double scale = 1.0;
    if      ( jetScale_ == +1 ) scale *= jet.shiftedEnUp();
    else if ( jetScale_ == -1 ) scale *= jet.shiftedEnDown();

    if ( isMC_ and !isSkipJER_ ) scale *= jet.smearedRes(jetResol_);

    return scale;
  }

  bool isGoodMuon(const cat::Muon& mu)
  {
    if ( std::abs(mu.eta()) > 2.1 ) return false;
    if ( std::isnan(mu.pt()) or mu.pt() < 27 ) return false;

    if ( mu.relIso(0.4) > 0.15 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const cat::Electron& el)
  {
    if ( std::abs(el.eta()) > 2.1 ) return false;
    if ( std::isnan(el.pt()) or el.pt() < 35 ) return false;

    if ( isMVAElectronSel_ and !el.isTrigMVAValid() ) return false;

    //if ( el.relIso(0.3) >= 0.11 ) return false;
    if ( !el.electronID(elIdName_) ) return false;
    //if ( !el.isPF() or !el.passConversionVeto() ) return false;
    const double scEta = std::abs(el.scEta());
    if ( isEcalCrackVeto_ and scEta > 1.4442 and scEta < 1.566 ) return false;
    //const double d0 = std::abs(el.dxy()), dz = std::abs(el.vz());
    //if      ( scEta <= 1.479 and (d0 > 0.05 or dz > 0.1) ) return false;
    //else if ( scEta >  1.479 and (d0 > 0.10 or dz > 0.2) ) return false;
    return true;
  }
  bool isVetoMuon(const cat::Muon& mu)
  {
    if ( std::abs(mu.eta()) > 2.4 ) return false;
    if ( std::isnan(mu.pt()) or mu.pt() < 10 ) return false;

    if ( !mu.isLooseMuon() ) return false;
    if ( mu.relIso(0.4) > 0.25 ) return false;
    return true;
  }
  bool isVetoElectron(const cat::Electron& el)
  {
    if ( std::abs(el.eta()) > 2.5 ) return false;
    if ( std::isnan(el.pt()) or el.pt() < 10 ) return false;
    if ( !el.electronID(elVetoIdName_) ) return false;
    //const double scEta = std::abs(el.scEta());
    //const double d0 = std::abs(el.dxy()), dz = std::abs(el.vz());
    //if      ( scEta <= 1.479 and (d0 > 0.05 or dz > 0.1) ) return false;
    //else if ( scEta >  1.479 and (d0 > 0.10 or dz > 0.2) ) return false;
    return true;
  }
  bool isBjet(const cat::Jet& jet)
  {
    using namespace cat;
    const double bTag = jet.bDiscriminator(bTagName_);
    if      ( bTagWP_ == BTagWP::CSVL ) return bTag > WP_BTAG_CSVv2L;
    else if ( bTagWP_ == BTagWP::CSVM ) return bTag > WP_BTAG_CSVv2M;
    else if ( bTagWP_ == BTagWP::CSVT ) return bTag > WP_BTAG_CSVv2T;
    return false;
  }

private:
  // Energy scales
  int muonScale_, electronScale_, jetScale_, jetResol_;
  bool isSkipJER_; // Do not apply JER, needed to remove randomness during the Synchronization

  // Efficiency SF
  cat::ScaleFactorEvaluator muonSF_, electronSF_;
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

private:
  TH1D* hWeight_;
  TH1D* hCutstep_, * hCutstepNoweight_;
  THnSparseF* hsp_;

  std::ofstream eventListFile_;

};

using namespace cat;

TopFCNCEventSelector::TopFCNCEventSelector(const edm::ParameterSet& pset):
  isMC_(pset.getParameter<bool>("isMC")),
  applyFilterAt_(pset.getParameter<int>("applyFilterAt"))
{
  const string eventFileName = pset.getUntrackedParameter<std::string>("eventFile", "");
  if ( ! eventFileName.empty() ) eventListFile_.open(eventFileName);

  const auto channel = pset.getParameter<std::string>("channel");
  if ( channel == "electron" ) channel_ = 11;
  else if ( channel == "muon" ) channel_ = 13;
  else edm::LogError("TopFCNCEventSelector") << "channel must be \"electron\" or \"muon\"\n";

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

  hCutstep_ = fs->make<TH1D>("hCutstep", "Cut step;;Events", nCutsteps, 0, nCutsteps);
  hCutstepNoweight_ = fs->make<TH1D>("hCutstepNoweight", "Cut step (no weight);;Events", nCutsteps, 0, nCutsteps);
  const char* labels[nCutsteps] = {
    "S1 all", "S2 goodVertex", "S3 trigger",
    "S4 lepton", "S5 other lepton veto", "S6 same lepton veto",
    "S7 nJet3"
  };
  for ( int i=0; i<nCutsteps; ++i ) {
    hCutstep_->GetXaxis()->SetBinLabel(i+1, labels[i]);
    hCutstepNoweight_->GetXaxis()->SetBinLabel(i+1, labels[i]);
  }

  constexpr int ndim = 1+1+2+4+2+2+6*5;
  const int nbins[ndim] = {nCutsteps, 50, 100, 100, // cutstep, st, ht
                           27, 100, 100, 100, // lepton pdgId, pt, eta, phi
                           100, 100, // met
                           10, 10, // nJet, nBjet
                           100, 100, 100, 25, 100,  // jet1
                           100, 100, 100, 25, 100,  // jet2
                           100, 100, 100, 25, 100,  // jet3
                           100, 100, 100, 25, 100,  // jet4
                           100, 100, 100, 25, 100,  // jet5
                           100, 100, 100, 25, 100}; // jet6
  const double pi = TMath::Pi();
  const double xmins[ndim] = {0, 0, 0, 0, -13.5, 0, -3, -pi, 0, -pi, 0, 0,
                              0, -3, -pi, 0, 0,  // jet1
                              0, -3, -pi, 0, 0,  // jet2
                              0, -3, -pi, 0, 0,  // jet3
                              0, -3, -pi, 0, 0,  // jet4
                              0, -3, -pi, 0, 0,  // jet5
                              0, -3, -pi, 0, 0}; // jet6
  const double xmaxs[ndim] = {nCutsteps, 50, 1000, 1000,
                              13.5, 500, 3, pi, 500, pi, 10, 10,
                              500, 3, pi, 50, 1,
                              500, 3, pi, 50, 1,
                              500, 3, pi, 50, 1,
                              500, 3, pi, 50, 1,
                              500, 3, pi, 50, 1,
                              500, 3, pi, 50, 1};
  const char* varNames[ndim] = {
    "cutstep", "vertex_n", "event_st", "event_ht", // [0-3]
    "lepton_pid", "lepton_pt", "lepton_eta", "lepton_phi", // [4-7]
    "met_pt", "met_phi", // [8,9]
    "jets_n", "bjets_n", // [10,11]
    "jet1_pt", "jet1_eta", "jet1_phi", "jet1_m", "jet1_btag", // [12-16]
    "jet2_pt", "jet2_eta", "jet2_phi", "jet2_m", "jet2_btag", // [17-21]
    "jet3_pt", "jet3_eta", "jet3_phi", "jet3_m", "jet3_btag", // [22-26]
    "jet4_pt", "jet4_eta", "jet4_phi", "jet4_m", "jet4_btag", // [27-31]
    "jet5_pt", "jet5_eta", "jet5_phi", "jet5_m", "jet5_btag", // [32-46]
    "jet6_pt", "jet6_eta", "jet6_phi", "jet6_m", "jet6_btag", // [37-41]
  };
  hsp_ = fs->make<THnSparseF>("hsp", "all infos", ndim, nbins, xmins, xmaxs);
  for ( int i=0; i<ndim; ++i ) {
    hsp_->GetAxis(i)->SetTitle(varNames[i]);
  }

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

  double vars[] = {0, 0, 0, 0,
    0, 0, 0, 0, // lepton
    0, 0, 0, 0, // met, jets_n, bjets_n
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // jet[1-6] pt, eta, phi, m, btag
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

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

  edm::Handle<int> nVertexHandle;
  event.getByToken(nVertexToken_, nVertexHandle);
  const int nVertex = *nVertexHandle;
  /*edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  const int nVertex = [&](){
    int n = 0;
    for ( auto& v : *vertexHandle ) {
      if ( v.isFake() or v.ndof() <= 4 ) continue;
      if ( std::abs(v.position().rho()) >= 2 or std::abs(v.z()) >= 24 ) continue;     
      ++n;
    }
    return n;
  }();*/
  vars[1] = nVertex;

  std::auto_ptr<std::vector<cat::Lepton> > out_leptons(new std::vector<cat::Lepton>());
  std::auto_ptr<std::vector<cat::Jet> > out_jets(new std::vector<cat::Jet>());

  // Compute event weight - from generator, pileup, etc
  double weight = 1.0;
  if ( isMC_ ) {
    edm::Handle<float> fHandle;
    edm::Handle<vfloat> vfHandle;

    float genWeight = 1.0;
    if ( genWeightIndex_ < 0 ) {
      event.getByToken(genWeightToken_, fHandle);
      genWeight = *fHandle;
    }
    else {
      event.getByToken(genWeightsToken_, vfHandle);
      genWeight = vfHandle->at(genWeightIndex_);
    }

    event.getByToken(pileupWeightToken_, fHandle);
    const float puWeight = *fHandle;

    weight *= genWeight*puWeight;
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

  double event_st = 0, event_ht = 0;

  // Select good leptons
  cat::MuonCollection selMuons, vetoMuons;
  for ( int i=0, n=muonHandle->size(); i<n; ++i ) {
    auto& p = muonHandle->at(i);
    const double scale = shiftedMuonScale(p);

    cat::Muon lep(p);
    lep.setP4(p.p4()*scale);
    if ( isGoodMuon(lep) ) selMuons.push_back(lep);
    else if ( isVetoMuon(lep) ) vetoMuons.push_back(lep);
    else continue;

    event_st += lep.pt();
    metDpx += lep.px()-p.px();
    metDpy += lep.py()-p.py();
  }
  cat::ElectronCollection selElectrons, vetoElectrons;
  for ( int i=0, n=electronHandle->size(); i<n; ++i ) {
    auto& p = electronHandle->at(i);
    const double scale = shiftedElectronScale(p)/p.smearedScale();

    cat::Electron lep(p);
    lep.setP4(p.p4()*scale);
    if ( isGoodElectron(lep) ) selElectrons.push_back(lep);
    else if ( isVetoElectron(lep) ) vetoElectrons.push_back(lep);
    else continue;

    event_st += lep.pt();
    metDpx += lep.px()-p.px()/p.smearedScale();
    metDpy += lep.py()-p.py()/p.smearedScale();
  }
  std::sort(selElectrons.begin(), selElectrons.end(),
            [&](const cat::Electron& a, const cat::Electron& b){ return a.pt() > b.pt(); });
  std::sort(selMuons.begin(), selMuons.end(),
            [&](const cat::Muon& a, const cat::Muon& b){ return a.pt() > b.pt(); });

  //const int leptons_n = selMuons.size() + selElectrons.size();
  const cat::Lepton* lepton1 = 0;
  double trigSF = 1, leptonSF = 1;
  if ( channel_ == 11 and !selElectrons.empty() ) {
    const auto& el = selElectrons.at(0);
    lepton1 = &el;
    leptonSF = electronSF_(el.pt(), std::abs(el.scEta()), electronSFShift_);
  }
  else if ( channel_ == 13 and !selMuons.empty() ) {
    const auto& mu = selMuons.at(0);
    lepton1 = &mu;
    leptonSF = muonSF_(mu.pt(), std::abs(mu.eta()), muonSFShift_);
  }
  if ( lepton1 ) {
    vars[4] = lepton1->pdgId();
    vars[5] = lepton1->pt();
    vars[6] = lepton1->eta();
    vars[7] = lepton1->phi();
    out_leptons->push_back(*lepton1);
  }

  // Select good jets
  int bjets_n = 0;
  for ( int i=0, n=jetHandle->size(); i<n; ++i ) {
    auto& p = jetHandle->at(i);
    if ( std::abs(p.eta()) > 2.4 ) continue;
    if ( !p.LooseId() ) continue;

    const double scale = shiftedJetScale(p);
    cat::Jet jet(p);
    jet.setP4(scale*p.p4());

    metDpx += jet.px()-p.px();
    metDpy += jet.py()-p.py();
    if ( jet.pt() < 30 ) continue;

    if ( lepton1 and deltaR(jet.p4(), lepton1->p4()) < 0.4 ) continue;

    event_ht += jet.pt();
    if ( isBjet(p) ) ++bjets_n;

    out_jets->push_back(jet);
  }
  event_st += event_ht;
  const int jets_n = out_jets->size();
  std::sort(out_jets->begin(), out_jets->end(),
            [&](const cat::Jet& a, const cat::Jet& b){return a.pt() > b.pt();});

  // Update & calculate met
  const double met_pt = hypot(metP4.px()-metDpx, metP4.py()-metDpy);
  const double met_phi = atan2(metP4.py()-metDpy, metP4.px()-metDpx);
  vars[8] = met_pt;
  vars[9] = met_phi;

  // Check cut steps
  std::vector<bool> cutsteps(nCutsteps);
  cutsteps[0] = true; // always true
  cutsteps[1] = (nVertex > 0);
  if ( channel_ == 11 ) {
    cutsteps[2] = (!isIgnoreTrig_ and isTrigEl != 0);
    cutsteps[3] = (selElectrons.size() == 1);
    cutsteps[4] = (selMuons.size()+vetoMuons.size() == 0);
    cutsteps[5] = vetoElectrons.empty();
  }
  else if ( channel_ == 13 ) {
    cutsteps[2] = (!isIgnoreTrig_ and isTrigMu != 0);
    cutsteps[3] = (selMuons.size() == 1);
    cutsteps[4] = (selElectrons.size()+vetoElectrons.size() == 0);
    cutsteps[5] = vetoMuons.empty();
  }
  cutsteps[6] = (jets_n >= 3);

  // Run though the cut steps
  int cutstep = 0;
  double w = weight;
  for ( cutstep = 0; cutstep < nCutsteps and cutsteps[cutstep]; ++cutstep ) {
    if ( cutstep == 2 and !isIgnoreTrig_ ) w *= trigSF;
    else if ( cutstep == 3 ) w *= leptonSF;
    hCutstep_->Fill(cutstep, w);
    hCutstepNoweight_->Fill(cutstep);
  }

  // Fill n-dim histogram
  vars[2] = event_st;
  vars[3] = event_ht;
  vars[0] = cutstep;
  vars[10] = jets_n;
  vars[11] = bjets_n;
  for ( unsigned int i=0, n=std::min(6, jets_n); i<n; ++i ) {
    const auto& jet = out_jets->at(i);
    const int ii = 12+i*5;
    vars[ii+0] = jet.pt();
    vars[ii+1] = jet.eta();
    vars[ii+2] = jet.phi();
    vars[ii+3] = jet.mass();
    vars[ii+4] = jet.bDiscriminator(bTagName_);
  }
  hsp_->Fill(vars, w);

  event.put(std::auto_ptr<int>(new int((int)channel_)), "channel");
  event.put(std::auto_ptr<float>(new float(weight)), "weight");
  event.put(std::auto_ptr<float>(new float(met_pt)), "met");
  event.put(std::auto_ptr<float>(new float(met_phi)), "metphi");
  event.put(out_leptons, "leptons");
  event.put(out_jets, "jets");

  // Apply filter at the given step.
  if ( cutstep >= applyFilterAt_ ) {
    if ( eventListFile_.is_open() ) {
      const int run = event.id().run();
      const int lum = event.id().luminosityBlock();
      const int evt = event.id().event();
      char buffer[101];
      snprintf(buffer, 100, "%6d %6d %10d", run, lum, evt);
      eventListFile_ << buffer;
      snprintf(buffer, 100, "  %+2d  %6.2f %+4.2f %+4.2f", int(vars[4]), vars[5], vars[6], vars[7]);
      eventListFile_ << buffer;
      snprintf(buffer, 100, "    %6.1f  %+4.2f", vars[8], vars[9]);
      eventListFile_ << buffer;
      snprintf(buffer, 100, "    %d %d", int(vars[10]), int(vars[11]));
      eventListFile_ << buffer;
      snprintf(buffer, 100, "  %6.2f %+4.2f %+4.2f  \n", vars[12], vars[13], vars[14]);
      eventListFile_ << buffer;
    }

    return true;
  }

  return false;
}

TopFCNCEventSelector::~TopFCNCEventSelector()
{
  if ( hCutstepNoweight_ ) {
    cout << "---- cut flows without weight ----\n";
    if ( channel_ == 11 ) cout << "Electron channel:\n";
    else if ( channel_ == 13 ) cout << "Muon channel:\n";
    size_t fw = 8;
    for ( int i=1; i<=nCutsteps; ++i ) {
      fw = std::max(fw, strlen(hCutstepNoweight_->GetXaxis()->GetBinLabel(i)));
    }
    fw += 2;
    for ( int i=1; i<=nCutsteps; ++i ) {
      const string name(hCutstepNoweight_->GetXaxis()->GetBinLabel(i));
      cout << name;
      for ( size_t k=name.size(); k<fw; ++k ) cout << ' ';
      cout << hCutstepNoweight_->GetBinContent(i) << '\n';
    }
    cout << "-----------------------------------\n";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopFCNCEventSelector);

