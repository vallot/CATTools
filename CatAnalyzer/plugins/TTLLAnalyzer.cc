#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "CATTools/DataFormats/interface/Lepton.h"
#include "CATTools/DataFormats/interface/Jet.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "CATTools/DataFormats/interface/GenTop.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"

#include "TTree.h"
#include "TString.h"

#include <memory>

class TTLLAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks>
{
public:
  TTLLAnalyzer(const edm::ParameterSet& pset);
  ~TTLLAnalyzer() {};

  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&);
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) {};
  void analyze(const edm::Event& event, const edm::EventSetup&) override;

private:
  bool isMC_, isTTbar_, doGenWeightSysts_;

  typedef std::vector<float> vfloat;
  edm::EDGetTokenT<float> topPtWeightToken_;
  edm::EDGetTokenT<float> puWeightToken_, puUpWeightToken_, puDnWeightToken_;
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<vfloat> scaleUpWeightsToken_, scaleDnWeightsToken_;
  edm::EDGetTokenT<vfloat> pdfWeightsToken_;
  edm::EDGetTokenT<float> csvWeightToken_;
  edm::EDGetTokenT<vfloat> csvWeightSystToken_;
  edm::EDGetTokenT<int> nVertexToken_;

  edm::EDGetTokenT<cat::LeptonCollection> leptonToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<float> metPtToken_, metPhiToken_;

  edm::EDGetTokenT<int> partonChToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTopToken_;
  edm::EDGetTokenT<std::vector<int>> partonModesToken_;

private:
  void clear();
  std::unique_ptr<cat::KinematicSolver> solver_;

  TTree* tree_;
  unsigned int b_run;
  unsigned short b_lumi;
  unsigned long long b_event;

  float b_wgt_gen;
  float b_wgt_pu, b_wgt_puUp, b_wgt_puDn;
  float b_wgt_csv;
  float b_wgt_topPt;

  std::unique_ptr<std::vector<float>> b_wgts_pdf;
  std::unique_ptr<std::vector<float>> b_wgts_scaleUp, b_wgts_scaleDn;
  std::unique_ptr<std::vector<float>> b_wgts_csvUp, b_wgts_csvDn;

  float b_event_st;
  unsigned char b_vertex_n;

  const unsigned static int kMaxNLeptons = 10;
  unsigned int b_leptons_n;
  float b_leptons_pt[kMaxNLeptons], b_leptons_eta[kMaxNLeptons], b_leptons_phi[kMaxNLeptons];
  float b_leptons_iso[kMaxNLeptons];
  unsigned char b_leptons_pid[kMaxNLeptons];

  float b_met_pt, b_met_phi;

  const unsigned static int kMaxNJets = 100;
  unsigned int b_jets_n;
  float b_jets_ht;
  float b_jets_pt[kMaxNJets], b_jets_eta[kMaxNJets], b_jets_phi[kMaxNJets];
  float b_jets_m[kMaxNJets];
  float b_jets_csv[kMaxNJets], b_jets_CvsL[kMaxNJets], b_jets_CvsB[kMaxNJets];
  unsigned char b_jets_hadFlav[kMaxNJets];

  unsigned char b_bjetsL_n, b_bjetsM_n, b_bjetsT_n;

  float b_ll_m, b_ll_pt, b_ll_eta, b_ll_phi;
  char b_ll_q;

  float b_kinTop1_m, b_kinTop1_pt, b_kinTop1_y;
  float b_kinTop2_m, b_kinTop2_pt, b_kinTop2_y;
  float b_kinTop1_csv, b_kinTop2_csv;
  float b_kinTT_m, b_kinTT_pt, b_kinTT_y, b_kinTT_dphi;

  const unsigned static int kMaxNPartons = 10;
  unsigned int b_tops_n;
  unsigned char b_tops_mode[kMaxNPartons];
  float b_tops_pt[kMaxNPartons], b_tops_eta[kMaxNPartons], b_tops_phi[kMaxNPartons], b_tops_m[kMaxNPartons];
  float b_ws_pt[kMaxNPartons], b_ws_eta[kMaxNPartons], b_ws_phi[kMaxNPartons], b_ws_m[kMaxNPartons];
  float b_wdaus1_pt[kMaxNPartons], b_wdaus1_eta[kMaxNPartons], b_wdaus1_phi[kMaxNPartons], b_wdaus1_m[kMaxNPartons];
  float b_wdaus2_pt[kMaxNPartons], b_wdaus2_eta[kMaxNPartons], b_wdaus2_phi[kMaxNPartons], b_wdaus2_m[kMaxNPartons];
  char b_wdaus1_id[kMaxNPartons], b_wdaus2_id[kMaxNPartons];
};

using namespace std;

TTLLAnalyzer::TTLLAnalyzer(const edm::ParameterSet& pset)
{
  isMC_ = pset.getParameter<bool>("isMC");
  isTTbar_ = pset.getParameter<bool>("isTTbar");
  doGenWeightSysts_ = pset.getParameter<bool>("doGenWeightSysts");

  auto srcTag = pset.getParameter<edm::InputTag>("src");
  leptonToken_ = consumes<cat::LeptonCollection>(edm::InputTag(srcTag.label(), "leptons"));
  jetToken_ = consumes<cat::JetCollection>(edm::InputTag(srcTag.label(), "jets"));
  metPtToken_  = consumes<float>(edm::InputTag(srcTag.label(), "met"));
  metPhiToken_ = consumes<float>(edm::InputTag(srcTag.label(), "metphi"));

  nVertexToken_ = consumes<int>(pset.getParameter<edm::InputTag>("nVertex"));
  if ( isMC_ ) {
    auto puWeightLabel = pset.getParameter<edm::InputTag>("puWeight");
    puWeightToken_ = consumes<float>(puWeightLabel);
    puUpWeightToken_ = consumes<float>(edm::InputTag(puWeightLabel.label(), "up"));
    puDnWeightToken_ = consumes<float>(edm::InputTag(puWeightLabel.label(), "dn"));

    auto genWeightLabel = pset.getParameter<edm::InputTag>("genWeight");
    genWeightToken_ = consumes<float>(genWeightLabel);
    scaleUpWeightsToken_ = consumes<vfloat>(edm::InputTag(genWeightLabel.label(), "scaleup"));
    scaleDnWeightsToken_ = consumes<vfloat>(edm::InputTag(genWeightLabel.label(), "scaledown"));
    pdfWeightsToken_ = consumes<vfloat>(edm::InputTag(genWeightLabel.label(), "pdf"));

    auto csvWeightLabel = pset.getParameter<edm::InputTag>("csvWeight");
    csvWeightToken_ = consumes<float>(edm::InputTag(csvWeightLabel.label()));
    csvWeightSystToken_ = consumes<vfloat>(edm::InputTag(csvWeightLabel.label(), "syst"));

    if ( isTTbar_ ) {
      topPtWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("topPtWeight"));

      auto partonLabel = pset.getParameter<edm::InputTag>("partonTop");
      partonChToken_ = consumes<int>(edm::InputTag(partonLabel.label(), "channel"));
      partonTopToken_ = consumes<reco::GenParticleCollection>(partonLabel);
      partonModesToken_ = consumes<std::vector<int>>(edm::InputTag(partonLabel.label(), "modes"));
    }
  }
  auto solverPSet = pset.getParameter<edm::ParameterSet>("solver");
  solver_.reset(new cat::DESYSmearedSolver(solverPSet));

  b_wgts_scaleUp.reset(new std::vector<float>());
  b_wgts_scaleDn.reset(new std::vector<float>());
  b_wgts_pdf.reset(new std::vector<float>());
  b_wgts_csvUp.reset(new std::vector<float>());
  b_wgts_csvDn.reset(new std::vector<float>());

  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("events", "TTLL event tree");

  tree_->Branch("run", &b_run, "run/i"); // run number in 10^7
  tree_->Branch("lumi", &b_lumi, "lumi/s"); // enough with 65535
  tree_->Branch("event", &b_event, "event/l"); // event number with 64bit int

  if ( isMC_ ) {
    // event weights
    tree_->Branch("wgt_gen", &b_wgt_gen, "wgt_gen/F"); // central gen weight
    tree_->Branch("wgt_pu", &b_wgt_pu, "wgt_pu/F"); // pileup weight
    tree_->Branch("wgt_csv", &b_wgt_csv, "wgt_csv/F"); // CSV weight
    if ( isTTbar_ ) {
      tree_->Branch("wgt_topPt", &b_wgt_topPt, "wgt_topPt/F");// top pt weight

      tree_->Branch("tops_n", &b_tops_n, "tops_n/i"); // array index must be an integer

      tree_->Branch("tops_mode", &b_tops_mode, "tops_mode[tops_n]/b");

      tree_->Branch("tops_pt" , &b_tops_pt , "tops_pt[tops_n]/F" );
      tree_->Branch("tops_eta", &b_tops_eta, "tops_eta[tops_n]/F");
      tree_->Branch("tops_phi", &b_tops_phi, "tops_phi[tops_n]/F");
      tree_->Branch("tops_m"  , &b_tops_m  , "tops_m[tops_n]/F"  );

      tree_->Branch("ws_pt" , &b_ws_pt , "ws_pt[tops_n]/F" );
      tree_->Branch("ws_eta", &b_ws_eta, "ws_eta[tops_n]/F");
      tree_->Branch("ws_phi", &b_ws_phi, "ws_phi[tops_n]/F");
      tree_->Branch("ws_m"  , &b_ws_m  , "ws_m[tops_n]/F"  );

      tree_->Branch("wdaus1_pt" , &b_wdaus1_pt , "wdaus1_pt[tops_n]/F" );
      tree_->Branch("wdaus1_eta", &b_wdaus1_eta, "wdaus1_eta[tops_n]/F");
      tree_->Branch("wdaus1_phi", &b_wdaus1_phi, "wdaus1_phi[tops_n]/F");
      tree_->Branch("wdaus1_m"  , &b_wdaus1_m  , "wdaus1_m[tops_n]/F"  );

      tree_->Branch("wdaus2_pt" , &b_wdaus2_pt , "wdaus2_pt[tops_n]/F" );
      tree_->Branch("wdaus2_eta", &b_wdaus2_eta, "wdaus2_eta[tops_n]/F");
      tree_->Branch("wdaus2_phi", &b_wdaus2_phi, "wdaus2_phi[tops_n]/F");
      tree_->Branch("wdaus2_m"  , &b_wdaus2_m  , "wdaus2_m[tops_n]/F"  );

      tree_->Branch("wdaus1_id", &b_wdaus1_id, "wdaus1_id[tops_n]/B");
      tree_->Branch("wdaus2_id", &b_wdaus2_id, "wdaus2_id[tops_n]/B");
    }

    tree_->Branch("wgt_puUp", &b_wgt_puUp, "wgt_puUp/F"); // pileup weight
    tree_->Branch("wgt_puDn", &b_wgt_puDn, "wgt_puDn/F"); // pileup weight

    if ( doGenWeightSysts_ ) {
      tree_->Branch("wgts_scaleUp", "std::vector<float>", &*b_wgts_scaleUp);
      tree_->Branch("wgts_scaleDn", "std::vector<float>", &*b_wgts_scaleDn);
      tree_->Branch("wgts_pdf", "std::vector<float>", &*b_wgts_pdf);
    }

    tree_->Branch("wgts_csvUp", "std::vector<float>", &*b_wgts_csvUp); // CSV weight
    tree_->Branch("wgts_csvDn", "std::vector<float>", &*b_wgts_csvDn); // CSV weight
  }

  // Event-wise observables
  tree_->Branch("event_st", &b_event_st, "event_st/F");
  tree_->Branch("vertex_n", &b_vertex_n, "vertex_n/b"); // enough with 255

  // Leptons
  tree_->Branch("leptons_n", &b_leptons_n, "leptons_n/i");
  tree_->Branch("leptons_pt" , &b_leptons_pt , "leptons_pt[leptons_n]/F");
  tree_->Branch("leptons_eta", &b_leptons_eta, "leptons_eta[leptons_n]/F");
  tree_->Branch("leptons_phi", &b_leptons_phi, "leptons_phi[leptons_n]/F");
  tree_->Branch("leptons_pid", &b_leptons_pid, "leptons_pid[leptons_n]/B"); // +-11 or +-13. enough with += 127
  tree_->Branch("leptons_iso", &b_leptons_iso, "leptons_iso[leptons_n]/F");

  // MET
  tree_->Branch("met_pt" , &b_met_pt , "met_pt/F");
  tree_->Branch("met_phi", &b_met_phi, "met_phi/F");

  // Jets
  tree_->Branch("jets_n", &b_jets_n, "jets_n/i"); // array index must be an integer
  tree_->Branch("jets_ht", &b_jets_ht, "jets_ht/F");
  tree_->Branch("jets_pt" , &b_jets_pt , "jets_pt[jets_n]/F");
  tree_->Branch("jets_eta", &b_jets_eta, "jets_eta[jets_n]/F");
  tree_->Branch("jets_phi", &b_jets_phi, "jets_phi[jets_n]/F");
  tree_->Branch("jets_m"  , &b_jets_m  , "jets_m[jets_n]/F");
  tree_->Branch("jets_csv" , &b_jets_csv , "jets_csv[jets_n]/F");
  tree_->Branch("jets_CvsL", &b_jets_CvsL, "jets_CvsL[jets_n]/F");
  tree_->Branch("jets_CvsB", &b_jets_CvsB, "jets_CvsB[jets_n]/F");
  if ( isMC_ ) tree_->Branch("jets_hadFlav", &b_jets_hadFlav, "jets_hadFlav[jets_n]/b");

  // bJets
  tree_->Branch("bjetsL_n", &b_bjetsL_n, "bjetsL_n/b"); // enough with 255
  tree_->Branch("bjetsM_n", &b_bjetsM_n, "bjetsM_n/b"); // enough with 255
  tree_->Branch("bjetsT_n", &b_bjetsT_n, "bjetsT_n/b"); // enough with 255

  // Composite objects
  tree_->Branch("ll_m"  , &b_ll_m  , "ll_m/F"  );
  tree_->Branch("ll_pt" , &b_ll_pt , "ll_pt/F" );
  tree_->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  tree_->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  tree_->Branch("ll_q"  , &b_ll_q  , "ll_q/B"  );

  tree_->Branch("kinTop1_m" , &b_kinTop1_m , "kinTop1_m/F" );
  tree_->Branch("kinTop1_pt", &b_kinTop1_pt, "kinTop1_pt/F");
  tree_->Branch("kinTop1_y" , &b_kinTop1_y , "kinTop1_y/F" );
  tree_->Branch("kinTop2_m" , &b_kinTop2_m , "kinTop2_m/F" );
  tree_->Branch("kinTop2_pt", &b_kinTop2_pt, "kinTop2_pt/F");
  tree_->Branch("kinTop2_y" , &b_kinTop2_y , "kinTop2_y/F" );
  tree_->Branch("kinTop1_csv", &b_kinTop1_csv, "kinTop1_csv/F");
  tree_->Branch("kinTop2_csv", &b_kinTop2_csv, "kinTop2_csv/F");

  tree_->Branch("kinTT_m"   , &b_kinTT_m   , "kinTT_m/F"   );
  tree_->Branch("kinTT_pt " , &b_kinTT_pt  , "kinTT_pt/F"  );
  tree_->Branch("kinTT_y"   , &b_kinTT_y   , "kinTT_y/F"   );
  tree_->Branch("kinTT_dphi", &b_kinTT_dphi, "kinTT_dphi/F");

}

void TTLLAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&)
{
  if ( dynamic_cast<cat::DESYSmearedSolver*>(solver_.get()) != 0 ) {
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(lumi.index());
    dynamic_cast<cat::DESYSmearedSolver*>(solver_.get())->setRandom(&engine);
  }
}

void TTLLAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&)
{
  edm::Handle<float> fHandle;
  edm::Handle<int> iHandle;
  edm::Handle<std::vector<float>> vfHandle;

  clear();
  b_run   = event.id().run();
  b_lumi  = event.id().luminosityBlock();
  b_event = event.id().event();

  if ( isMC_ ) {
    if ( isTTbar_ and event.getByToken(topPtWeightToken_, fHandle) ) {
      b_wgt_topPt = *fHandle;
    }

    event.getByToken(puWeightToken_, fHandle);
    b_wgt_pu = *fHandle;

    event.getByToken(puUpWeightToken_, fHandle);
    b_wgt_puUp = *fHandle;

    event.getByToken(puDnWeightToken_, fHandle);
    b_wgt_puDn = *fHandle;

    event.getByToken(genWeightToken_, fHandle);
    b_wgt_gen = *fHandle;

    if ( doGenWeightSysts_ and event.getByToken(scaleUpWeightsToken_, vfHandle) ) {
      for ( auto x : *vfHandle ) b_wgts_scaleUp->push_back(x);
    }
    if ( doGenWeightSysts_ and event.getByToken(scaleDnWeightsToken_, vfHandle) ) {
      for ( auto x : *vfHandle ) b_wgts_scaleDn->push_back(x);
    }

    if ( doGenWeightSysts_ and event.getByToken(pdfWeightsToken_, vfHandle) ) {
      for ( auto x : *vfHandle ) b_wgts_pdf->push_back(x);
    }

    event.getByToken(csvWeightToken_, fHandle);
    b_wgt_csv = *fHandle;
    event.getByToken(csvWeightSystToken_, vfHandle);
    for ( int i=0, n=vfHandle->size(); i<n; i+=2 ) {
      b_wgts_csvUp->push_back(vfHandle->at(i+0));
      b_wgts_csvDn->push_back(vfHandle->at(i+1));
    }
  }

  event.getByToken(nVertexToken_, iHandle);
  b_vertex_n = *iHandle;

  event.getByToken(metPtToken_, fHandle);
  b_met_pt  = *fHandle;
  b_event_st += b_met_pt;

  event.getByToken(metPhiToken_, fHandle);
  b_met_phi = *fHandle;

  edm::Handle<cat::LeptonCollection> leptonHandle;
  event.getByToken(leptonToken_, leptonHandle);
  b_leptons_n = leptonHandle->size();
  for ( unsigned int i=0; i<b_leptons_n; ++i ) {
    const auto& lepton = leptonHandle->at(i);
    if ( i < kMaxNLeptons ) {
      b_leptons_pt[i] = lepton.pt();
      b_leptons_eta[i] = lepton.eta();
      b_leptons_phi[i] = lepton.phi();
      b_leptons_pid[i] = lepton.pdgId();
      b_leptons_iso[i] = lepton.relIso(abs(lepton.pdgId()) == 13 ? 0.4 : 0.3);
    }

    b_event_st += lepton.pt();
  }

  edm::Handle<cat::JetCollection> jetHandle;
  event.getByToken(jetToken_, jetHandle);
  b_jets_n = jetHandle->size();
  for ( unsigned int i=0; i<b_jets_n; ++i ) {
    const auto& jet = jetHandle->at(i);
    if ( i < kMaxNJets ) {
      b_jets_pt[i] = jet.pt();
      b_jets_eta[i] = jet.eta();
      b_jets_phi[i] = jet.phi();
      b_jets_m[i] = jet.mass();
      b_jets_csv[i] = jet.bDiscriminator(cat::BTAG_CSVv2);
      b_jets_CvsL[i] = jet.bDiscriminator(cat::CTAG_CvsL);
      b_jets_CvsB[i] = jet.bDiscriminator(cat::CTAG_CvsB);
      if ( isMC_ ) b_jets_hadFlav[i] = jet.hadronFlavour();
    }
    b_jets_ht += jet.pt();

    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2L ) ++b_bjetsL_n;
    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2M ) ++b_bjetsM_n;
    if ( jet.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2T ) ++b_bjetsT_n;
  }

  switch ( 1 ) default: {
    if ( b_leptons_n < 2 ) break;
    const auto lep1P4 = leptonHandle->at(0).p4();
    const auto lep2P4 = leptonHandle->at(1).p4();
    const auto llP4 = lep1P4+lep2P4;
    b_ll_m = llP4.mass();
    b_ll_pt = llP4.pt();
    b_ll_eta = llP4.eta();
    b_ll_phi = llP4.phi();
    b_ll_q = leptonHandle->at(0).charge()+leptonHandle->at(1).charge();

    if ( b_jets_n < 2 ) break;

    const math::XYZTLorentzVector metP4(b_met_pt*cos(b_met_phi), b_met_pt*sin(b_met_phi), 0, b_met_pt);
    cat::KinematicSolution solution;
    double quality = -1;
    for ( unsigned int i=0; i<b_jets_n; ++i ) {
      auto& jet1 = jetHandle->at(i);
      const bool isBjet1 = ( jet1.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2M );
      for ( unsigned int j=0; j<b_jets_n; ++j ) {
        auto& jet2 = jetHandle->at(j);
        const bool isBjet2 = ( jet2.bDiscriminator(cat::BTAG_CSVv2) > cat::WP_BTAG_CSVv2M );

        if ( ( b_bjetsM_n >= 2 and (isBjet1 and isBjet2) ) or // use both b-jets for b 2 btag event
             ( b_bjetsM_n == 1 and (isBjet1 or isBjet2) ) or // use at least one b-jet for only one b tag event
             ( b_bjetsM_n == 0 ) ) { // doesnt' matter for 0 b tagged event
          solver_->solve(metP4, lep1P4, lep2P4, jet1.p4(), jet2.p4());
          if ( quality > solver_->quality() ) {
            quality = solver_->quality();
            solution = solver_->solution();
          }
        }
      }
    }
    if ( quality >= 0 ) {
      b_kinTop1_m  = solution.t1().mass();
      b_kinTop1_pt = solution.t1().pt();
      b_kinTop1_y  = solution.t1().Rapidity();

      b_kinTop2_m  = solution.t2().mass();
      b_kinTop2_pt = solution.t2().pt();
      b_kinTop2_y  = solution.t2().Rapidity();

      double minDR1 = 1e9, minDR2 = 1e9;
      unsigned int i1 = b_jets_n, i2 = b_jets_n;
      for ( unsigned int i=0; i<b_jets_n; ++i ) {
        const double dR1 = deltaR2(solution.j1().eta(), solution.j1().phi(), b_jets_eta[i], b_jets_phi[i]);
        const double dR2 = deltaR2(solution.j2().eta(), solution.j2().phi(), b_jets_eta[i], b_jets_phi[i]);
        if ( dR1 < minDR1 ) { minDR1 = dR1; i1 = i; }
        if ( dR2 < minDR2 ) { minDR2 = dR2; i2 = i; }
      }
      if ( i1 < b_jets_n ) b_kinTop1_csv = b_jets_csv[i1];
      if ( i2 < b_jets_n ) b_kinTop2_csv = b_jets_csv[i2];

      b_kinTT_m    = solution.tt().mass();
      b_kinTT_pt   = solution.tt().pt();
      b_kinTT_y    = solution.tt().Rapidity();
      b_kinTT_dphi = deltaPhi(solution.t1().phi(), solution.t2().phi());
    }
    
  }
  b_event_st += b_jets_ht;

  if ( isMC_ and isTTbar_ ) {
    edm::Handle<int> channelHandle;
    event.getByToken(partonChToken_, channelHandle);

    edm::Handle<std::vector<int>> modesHandle;
    event.getByToken(partonModesToken_, modesHandle);

    edm::Handle<reco::GenParticleCollection> partonTopHandle;
    event.getByToken(partonTopToken_, partonTopHandle);

    b_tops_n = modesHandle->size();
    for ( unsigned int i=0, n=std::min(b_tops_n, kMaxNPartons); i<n; ++i ) {
      b_tops_mode[i] = modesHandle->at(i);

      auto& top = partonTopHandle->at(i);
      b_tops_pt[i] = top.pt();
      b_tops_eta[i] = top.eta();
      b_tops_phi[i] = top.phi();
      b_tops_m[i] = top.mass();

      if ( top.numberOfDaughters() < 2 ) continue;
      auto& wboson = *top.daughter(0); // dau0 is W
      //auto& b1 = *top.daughter(1); // dau1 is b quark
      b_ws_pt[i] = wboson.pt();
      b_ws_eta[i] = wboson.eta();
      b_ws_phi[i] = wboson.phi();
      b_ws_m[i] = wboson.mass();

      if ( wboson.numberOfDaughters() < 2 ) continue;
      auto& wdau1 = *wboson.daughter(0); // Note: abs(dau1_id) < abs(dau2_id)
      auto& wdau2 = *wboson.daughter(1); // therefore, lepton comes first, neutrino last

      b_wdaus1_pt[i] = wdau1.pt();
      b_wdaus1_eta[i] = wdau1.eta();
      b_wdaus1_phi[i] = wdau1.phi();
      b_wdaus1_m[i] = wdau1.mass();
      b_wdaus1_id[i] = wdau1.pdgId();

      b_wdaus2_pt[i] = wdau2.pt();
      b_wdaus2_eta[i] = wdau2.eta();
      b_wdaus2_phi[i] = wdau2.phi();
      b_wdaus2_m[i] = wdau2.mass();
      b_wdaus2_id[i] = wdau2.pdgId();
    }
  }

  tree_->Fill();
}

void TTLLAnalyzer::clear()
{
  b_run = b_lumi = b_event = 0;

  b_wgt_gen = 1.;
  b_wgt_pu = b_wgt_puUp = b_wgt_puDn = 1.;
  b_wgt_csv = 1.;
  b_wgt_topPt = 1.;

  b_wgts_pdf->clear();
  b_wgts_scaleUp->clear();
  b_wgts_scaleDn->clear();
  b_wgts_csvUp->clear();
  b_wgts_csvDn->clear();

  b_event_st = -10;
  b_vertex_n = 0;
  b_leptons_n = 0;
  for ( unsigned int i=0; i<kMaxNLeptons; ++i ) {
    b_leptons_pt[i] = b_leptons_eta[i] = b_leptons_phi[i] = -10;
    b_leptons_pid[i] = 0;
  }
  b_met_pt = b_met_phi = -10;
  b_jets_n = 0;
  b_jets_ht = -10;
  for ( unsigned int i=0; i<kMaxNJets; ++i ) {
    b_jets_pt[i] = b_jets_eta[i] = b_jets_phi[i] = -10;
    b_jets_m[i] = -10;
    b_jets_csv[i] = b_jets_CvsL[i] = b_jets_CvsB[i] = -10;
    b_jets_hadFlav[i] = 0;
  }
  b_bjetsL_n = b_bjetsM_n = b_bjetsT_n = 0;

  b_ll_m = b_ll_pt = b_ll_eta = b_ll_phi = -10;
  b_ll_q = 0;
  b_kinTop1_m = b_kinTop1_pt = b_kinTop1_y = -10;
  b_kinTop2_m = b_kinTop2_pt = b_kinTop2_y = -10;
  b_kinTop1_csv = b_kinTop2_csv = -10;
  b_kinTT_m = b_kinTT_pt = b_kinTT_y = b_kinTT_dphi = -10;

  b_tops_n = 0;
  for ( unsigned int i=0; i<kMaxNPartons; ++i ) {
    b_tops_pt[i] = b_tops_eta[i] = b_tops_phi[i] = b_tops_m[i] = -10;
    b_ws_pt[i] = b_ws_eta[i] = b_ws_phi[i] = b_ws_m[i] = -10;
    b_wdaus1_pt[i] = b_wdaus1_eta[i] = b_wdaus1_phi[i] = b_wdaus1_m[i] = -10;
    b_wdaus2_pt[i] = b_wdaus2_eta[i] = b_wdaus2_phi[i] = b_wdaus2_m[i] = -10;
    b_wdaus1_id[i] = b_wdaus2_id[i] = 0;
  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLLAnalyzer);
