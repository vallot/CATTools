#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;
using namespace cat;

class TtbarDiLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer();

  enum {
    CH_NONE=0, CH_MUEL=1, CH_ELEL=2, CH_MUMU=3
  };

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  void selectMuons(const cat::MuonCollection& muons, ParticleCollection& selmuons) const;
  void selectElecs(const cat::ElectronCollection& elecs, ParticleCollection& selelecs) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const ParticleCollection& recolep) const;
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  const reco::Candidate* getLast(const reco::Candidate* p) const;

  edm::EDGetTokenT<int> recoFiltersToken_;
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<vector<int> > partonTop_modes_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTop_genParticles_;

  edm::EDGetTokenT<reco::GenParticleCollection> pseudoTop_;

  TTree * ttree_;
  int b_partonChannel, b_partonMode1, b_partonMode2;
  float b_partonlep1_pt, b_partonlep1_eta;
  float b_partonlep2_pt, b_partonlep2_eta;
  int b_pseudoTopChannel;
  float b_pseudoToplep1_pt, b_pseudoToplep1_eta;
  float b_pseudoToplep2_pt, b_pseudoToplep2_eta;
  int b_njet, b_nbjet, b_step, b_channel;
  bool b_lepinPhase, b_jetinPhase;
  float b_MET, b_maxweight;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;
  float b_jet1_pt, b_jet1_eta, b_jet1_CSVInclV2;
  float b_jet2_pt, b_jet2_eta, b_jet2_CSVInclV2;
  float b_top1_pt, b_top1_eta, b_top1_phi, b_top1_rapi;
  float b_top2_pt, b_top2_eta, b_top2_phi, b_top2_rapi;
  float b_ttbar_pt, b_ttbar_eta, b_ttbar_phi, b_ttbar_m, b_ttbar_rapi;
  int b_tri;
  int b_filtered;
  int b_is3lep;

  std::unique_ptr<TtFullLepKinSolver> solver;
  bool isTTbarMC_;
  //enum TTbarMode { CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON };
  //enum DecayMode { CH_HADRON = 1, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

  const static int NCutflow = 6;
  std::vector<std::vector<int> > cutflow_;
};
//
// constructors and destructor
//
TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  trigTokenMUEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));

  muonToken_ = consumes<cat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  isTTbarMC_ = iConfig.getParameter<bool>("isTTbarMC");
  if ( isTTbarMC_ ) {
    partonTop_channel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("partonTop_channel"));
    partonTop_modes_   = consumes<vector<int> >(iConfig.getParameter<edm::InputTag>("partonTop_modes"));
    partonTop_genParticles_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("partonTop_genParticles"));

    pseudoTop_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pseudoTop"));
  }

  const double tmassbegin = iConfig.getParameter<double>       ("tmassbegin");
  const double tmassend   = iConfig.getParameter<double>       ("tmassend");
  const double tmassstep  = iConfig.getParameter<double>       ("tmassstep");
  const auto   nupars     = iConfig.getParameter<vector<double> >("neutrino_parameters");

  solver.reset(new TtFullLepKinSolver(tmassbegin, tmassend, tmassstep, nupars));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("tree", "tree");
  ttree_->Branch("parton_channel", &b_partonChannel, "parton_channel/I");
  ttree_->Branch("parton_mode1", &b_partonMode1, "parton_mode1/I");
  ttree_->Branch("partonlep1_pt", &b_partonlep1_pt, "partonlep1_pt/F");
  ttree_->Branch("partonlep1_eta", &b_partonlep1_eta, "partonlep1_eta/F");
  ttree_->Branch("partonlep2_pt", &b_partonlep2_pt, "partonlep2_pt/F");
  ttree_->Branch("partonlep2_eta", &b_partonlep2_eta, "partonlep2_eta/F");
  ttree_->Branch("parton_mode2", &b_partonMode2, "parton_mode2/I");

  ttree_->Branch("pseudoTop_channel", &b_pseudoTopChannel, "pseudoTop_channel/I");
  ttree_->Branch("pseudoToplep1_pt", &b_pseudoToplep1_pt, "pseudoToplep1_pt/F");
  ttree_->Branch("pseudoToplep1_eta", &b_pseudoToplep1_eta, "pseudoToplep1_eta/F");
  ttree_->Branch("pseudoToplep2_pt", &b_pseudoToplep2_pt, "pseudoToplep2_pt/F");
  ttree_->Branch("pseudoToplep2_eta", &b_pseudoToplep2_eta, "pseudoToplep2_eta/F");

  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("nbjet", &b_nbjet, "nbjet/I");
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("step", &b_step, "step/I");
  ttree_->Branch("lepinPhase", &b_lepinPhase, "lepinPhase/O");
  ttree_->Branch("jetinPhase", &b_jetinPhase, "jetinPhase/O");

  ttree_->Branch("lep1_pt", &b_lep1_pt, "lep1_pt/F");
  ttree_->Branch("lep1_eta", &b_lep1_eta, "lep1_eta/F");
  ttree_->Branch("lep1_phi", &b_lep1_phi, "lep1_phi/F");
  ttree_->Branch("lep2_pt", &b_lep2_pt, "lep2_pt/F");
  ttree_->Branch("lep2_eta", &b_lep2_eta, "lep2_eta/F");
  ttree_->Branch("lep2_phi", &b_lep2_phi, "lep2_phi/F");
  ttree_->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
  ttree_->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  ttree_->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  ttree_->Branch("ll_m", &b_ll_m, "ll_m/F");
  ttree_->Branch("jet1_pt", &b_jet1_pt, "jet1_pt/F");
  ttree_->Branch("jet2_pt", &b_jet2_pt, "jet2_pt/F");
  ttree_->Branch("jet1_eta", &b_jet1_eta, "jet1_eta/F");
  ttree_->Branch("jet2_eta", &b_jet2_eta, "jet2_eta/F");
  ttree_->Branch("jet1_CSVInclV2", &b_jet1_CSVInclV2, "jet1_CSVInclV2/F");
  ttree_->Branch("jet2_CSVInclV2", &b_jet2_CSVInclV2, "jet2_CSVInclV2/F");

  ttree_->Branch("top1_pt", &b_top1_pt, "top1_pt/F");
  ttree_->Branch("top1_eta", &b_top1_eta, "top1_eta/F");
  ttree_->Branch("top1_phi", &b_top1_phi, "top1_phi/F");
  ttree_->Branch("top1_rapi", &b_top1_rapi, "top1_rapi/F");
  ttree_->Branch("top2_pt", &b_top2_pt, "top2_pt/F");
  ttree_->Branch("top2_eta", &b_top2_eta, "top2_eta/F");
  ttree_->Branch("top2_phi", &b_top2_phi, "top2_phi/F");
  ttree_->Branch("top2_rapi", &b_top2_rapi, "top2_rapi/F");
  ttree_->Branch("ttbar_pt", &b_ttbar_pt, "ttbar_pt/F");
  ttree_->Branch("ttbar_eta", &b_ttbar_eta, "ttbar_eta/F");
  ttree_->Branch("ttbar_phi", &b_ttbar_phi, "ttbar_phi/F");
  ttree_->Branch("ttbar_rapi", &b_ttbar_rapi, "ttbar_rapi/F");
  ttree_->Branch("ttbar_m", &b_ttbar_m, "ttbar_m/F");

  ttree_->Branch("tri", &b_tri, "tri/I");
  ttree_->Branch("filtered", &b_filtered, "filtered/I");
  ttree_->Branch("is3lep", &b_is3lep, "is3lep/I");

  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});
}

TtbarDiLeptonAnalyzer::~TtbarDiLeptonAnalyzer()
{
  cout <<"cut flow         emu         ee         mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step "<< i << " "<< cutflow_[i][0] <<  " "<< cutflow_[i][1] << " " << cutflow_[i][2] << " " << cutflow_[i][3]<< endl;
  }
}

void TtbarDiLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_partonChannel = -1; b_partonMode1 = -1; b_partonMode2 = -1;
  b_partonlep1_pt = -9; b_partonlep1_eta = -9;
  b_partonlep2_pt = -9; b_partonlep2_eta = -9;
  b_pseudoTopChannel = -1;
  b_pseudoToplep1_pt = -9; b_pseudoToplep1_eta = -9;
  b_pseudoToplep2_pt = -9; b_pseudoToplep2_eta = -9;
  b_MET = -1;
  b_njet = -1;
  b_nbjet = -1;
  b_channel = 0;
  b_step = -1;
  b_lepinPhase = false; b_jetinPhase = false;
  b_lep1_pt = -9; b_lep1_eta = -9; b_lep1_phi = -9;
  b_lep2_pt = -9; b_lep2_eta = -9; b_lep2_phi = -9;
  b_ll_pt = -9; b_ll_eta = -9; b_ll_phi = -9; b_ll_m = -9;
  b_jet1_pt = -9; b_jet1_eta = -9; b_jet1_CSVInclV2 = -9;
  b_jet2_pt = -9; b_jet2_eta = -9; b_jet2_CSVInclV2 = -9;
  b_top1_pt = -9; b_top1_eta = -9; b_top1_phi = -9; b_top1_rapi = -9;
  b_top2_pt = -9; b_top2_eta = -9; b_top2_phi = -9; b_top2_rapi = -9;
  b_ttbar_pt = -9; b_ttbar_eta = -9; b_ttbar_phi = -9; b_ttbar_m = -9; b_ttbar_rapi = -9;
  b_tri = -9;
  b_filtered = -9; b_is3lep = -9;
  if ( isTTbarMC_ and iEvent.isRealData() ) isTTbarMC_ = false;

  // bool debug = false;
  // if (iEvent.id().event() == 312909020 || iEvent.id().event() == 255013550){
  //   debug = true;
  //   cout <<"############## debugging iEvent.id().event()" << iEvent.id().event()<<endl;
  // }
  edm::Handle<reco::VertexCollection> vertices;      iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ return;} // skip the event if no PV found
  // const reco::Vertex &PV = vertices->front();
  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);

  if (isTTbarMC_){
    edm::Handle<int> partonTop_channel;
    edm::Handle<vector<int> > partonTop_modes;
    edm::Handle<reco::GenParticleCollection> partonTop_genParticles;
    iEvent.getByToken(partonTop_channel_, partonTop_channel);
    iEvent.getByToken(partonTop_modes_, partonTop_modes);
    iEvent.getByToken(partonTop_genParticles_, partonTop_genParticles);
    if ( (*partonTop_modes).size() == 0 ) {
      b_partonMode1 = 0;
      b_partonMode2 = 0;
    }
    else if ( (*partonTop_modes).size() == 1 ) { b_partonMode2 = 0; }
    else{
      b_partonChannel = *partonTop_channel;
      b_partonMode1 = (*partonTop_modes)[0];
      b_partonMode2 = (*partonTop_modes)[1];
    }

    if ( !(partonTop_genParticles->empty()) ){

      // Get Top quark pairs
      const auto parton1 = &partonTop_genParticles->at(0);
      const auto parton2 = &partonTop_genParticles->at(1);
      // Get W and b quarks
      if ( parton1 and parton2 ) {
        const auto partonW1 = parton1->daughter(0);
        const auto partonB1 = parton1->daughter(1);
        const auto partonW2 = parton2->daughter(0);
        const auto partonB2 = parton2->daughter(1);

        // Get W daughters
        if ( partonW1 and partonW2 and partonB1 and partonB2 ) {
          const auto partonW11 = partonW1->daughter(0);
          const auto partonW21 = partonW2->daughter(0);

          // Fill lepton informations
          b_partonlep1_pt = partonW11->pt();
          b_partonlep1_eta = partonW11->eta();
          b_partonlep2_pt = partonW21->pt();
          b_partonlep2_eta = partonW21->eta();
        }
      }
    }

    edm::Handle<reco::GenParticleCollection> pseudoTopHandle;
    iEvent.getByToken(pseudoTop_          , pseudoTopHandle);
    if ( !(pseudoTopHandle->empty()) ){
      b_pseudoTopChannel = CH_NONE;

      // Get Top quark pairs
      const auto pseudoTop1 = &pseudoTopHandle->at(0);
      const auto pseudoTop2 = &pseudoTopHandle->at(1);

      // Get W and b quarks
      if ( pseudoTop1 and pseudoTop2 ) {
        const auto pseudoW1 = pseudoTop1->daughter(0);
        const auto pseudoB1 = pseudoTop1->daughter(1);
        const auto pseudoW2 = pseudoTop2->daughter(0);
        const auto pseudoB2 = pseudoTop2->daughter(1);

        // Get W daughters
        if ( pseudoW1 and pseudoW2 and pseudoB1 and pseudoB2 ) {
          const auto pseudoW11 = pseudoW1->daughter(0);
          const auto pseudoW21 = pseudoW2->daughter(0);

          // Fill leps informations
          const int pseudoW1DauId = abs(pseudoW11->pdgId());
          const int pseudoW2DauId = abs(pseudoW21->pdgId());
          b_pseudoToplep1_pt = pseudoW11->pt();
          b_pseudoToplep1_eta = pseudoW11->eta();
          b_pseudoToplep2_pt = pseudoW21->pt();
          b_pseudoToplep2_eta = pseudoW21->eta();
          if ( pseudoW1DauId > 10 and pseudoW2DauId > 10 ) {
            switch ( pseudoW1DauId+pseudoW2DauId ) {
              case 22: b_pseudoTopChannel = CH_ELEL; break;
              case 26: b_pseudoTopChannel = CH_MUMU; break;
              default: b_pseudoTopChannel = CH_MUEL;
            }
          }
        }
      }
    }
  }

  // Store reco filter results
  edm::Handle<int> recoFiltersHandle;
  iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;

  // Find leptons and sort by pT
  ParticleCollection recolep;
  selectMuons(*muons, recolep);
  selectElecs(*electrons, recolep);
  if (recolep.size() < 2){
    ttree_->Fill();
    return;
  }
  sort(recolep.begin(), recolep.end(), GtByCandPt());
  const cat::Particle& recolep1 = recolep[0];
  const cat::Particle& recolep2 = recolep[1];

  // Determine channel
  const int pdgIdSum = std::abs(recolep1.pdgId()) + std::abs(recolep2.pdgId());
  if (pdgIdSum == 24) b_channel = CH_MUEL; // emu
  if (pdgIdSum == 22) b_channel = CH_ELEL; // ee
  if (pdgIdSum == 26) b_channel = CH_MUMU; // mumu

  // Trigger results
  edm::Handle<int> trigHandle;
  if      ( b_channel == CH_ELEL ) iEvent.getByToken(trigTokenELEL_, trigHandle);
  else if ( b_channel == CH_MUMU ) iEvent.getByToken(trigTokenMUMU_, trigHandle);
  else if ( b_channel == CH_MUEL ) iEvent.getByToken(trigTokenMUEL_, trigHandle);
  b_tri = *trigHandle;

  cutflow_[++b_step][b_channel]++;

  b_lep1_pt = recolep1.pt(); b_lep1_eta = recolep1.eta(); b_lep1_phi = recolep1.phi();
  b_lep2_pt = recolep2.pt(); b_lep2_eta = recolep2.eta(); b_lep2_phi = recolep2.phi();
  const TLorentzVector tlv_ll = recolep1.tlv()+recolep2.tlv();
  b_ll_pt = tlv_ll.Pt(); b_ll_eta = tlv_ll.Eta(); b_ll_phi = tlv_ll.Phi(); b_ll_m = tlv_ll.M();

  if (b_ll_m < 20.){
    ttree_->Fill();
    return;
  }
  if (recolep1.charge() * recolep2.charge() > 0){
    ttree_->Fill();
    return;
  }
  cutflow_[++b_step][b_channel]++;

  if (b_channel != CH_MUEL){
    if ((b_ll_m > 76) && (b_ll_m < 106)){
      ttree_->Fill();
      return;
    }
  }
  cutflow_[++b_step][b_channel]++;

  JetCollection&& selectedJets = selectJets(*jets, recolep);
  JetCollection&& selectedBJets = selectBJets(selectedJets);
  const TLorentzVector met = mets->front().tlv();
  b_MET = met.Pt();
  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();

  if (selectedJets.size() < 2){
    ttree_->Fill();
    return;
  }
  cutflow_[++b_step][b_channel]++;

  if (b_channel != CH_MUEL){
    if (b_MET < 40.){
      ttree_->Fill();
      return;
    }
  }
  cutflow_[++b_step][b_channel]++;

  if (selectedBJets.size() == 0){
    ttree_->Fill();
    return;
  }
  cutflow_[++b_step][b_channel]++;

  ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
  //int kin=0;
  TLorentzVector top1, top2, nu1, nu2;
  double maxweight=0;
  //const cat::Jet* kinj1, * kinj2;

  const TLorentzVector recolepLV1= recolep1.tlv();
  const TLorentzVector recolepLV2= recolep2.tlv();
  for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
    const TLorentzVector recojet1= jet1->tlv();
    for (auto jet2 = next(jet1); jet2 != end; ++jet2){

      const TLorentzVector recojet2= jet2->tlv();

      b_jet1_pt = recojet1.Pt();
      b_jet1_eta = recojet1.Eta();
      b_jet2_pt = recojet2.Pt();
      b_jet2_eta = recojet2.Eta();

      const double xconstraint = recolep1.px()+recolep2.px()+ recojet1.Px() + recojet2.Px() +met.Px();
      const double yconstraint = recolep1.py()+recolep2.py()+ recojet1.Py() + recojet2.Py() +met.Py();

      solver->SetConstraints(xconstraint, yconstraint);
      const auto nuSol1 = solver->getNuSolution( recolepLV1, recolepLV2 , recojet1, recojet2);
      const auto nuSol2 = solver->getNuSolution( recolepLV1, recolepLV2 , recojet2, recojet1);

      const double weight1 = nuSol1.weight;
      const double weight2 = nuSol2.weight;

      if ( weight1 > maxweight and weight1 >= weight2 ) {
        nu1 = cat::ToTLorentzVector(nuSol1.neutrino);
        nu2 = cat::ToTLorentzVector(nuSol1.neutrinoBar);
        maxweight = weight1;
      }
      else if ( weight2 > maxweight and weight2 >= weight1 ) {
        nu1 = cat::ToTLorentzVector(nuSol2.neutrino);
        nu2 = cat::ToTLorentzVector(nuSol2.neutrinoBar);
        maxweight = weight2;
      }
      else continue;

      top1 = recolepLV1+recojet1+nu1;
      top2 = recolepLV2+recojet2+nu2;
    }
  }

  b_top1_pt = top1.Pt();
  b_top1_eta = top1.Eta();
  b_top1_phi = top1.Phi();
  b_top1_rapi = top1.Rapidity();
  b_top2_pt = top2.Pt();
  b_top2_eta = top2.Eta();
  b_top2_phi = top2.Phi();
  b_top2_rapi = top2.Rapidity();

  TLorentzVector ttbar = top1+top2;
  b_ttbar_pt = ttbar.Pt();
  b_ttbar_eta = ttbar.Eta();
  b_ttbar_phi = ttbar.Phi();
  b_ttbar_m = ttbar.M();
  b_ttbar_rapi = ttbar.Rapidity();

  b_maxweight = maxweight;
  //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
  // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_MET, b_ll_mass, b_maxweight);

  ttree_->Fill();
}

const reco::Candidate* TtbarDiLeptonAnalyzer::getLast(const reco::Candidate* p) const
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

void TtbarDiLeptonAnalyzer::selectMuons(const cat::MuonCollection& muons, ParticleCollection& selmuons) const
{
  for (auto& mu : muons) {
    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    //if (!mu.isMediumMuon()) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.12) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }
}

void TtbarDiLeptonAnalyzer::selectElecs(const cat::ElectronCollection& elecs, ParticleCollection& selelecs) const
{
  for (auto& el : elecs) {
    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    //if (!el.electronID("cutBasedElectronID-Spring15-50ns-V1-standalone-medium")) continue;
    //if (el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium") == 0) continue;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) continue;
    //if (!el.passConversionVeto()) continue;
    //if (!el.isPF()) continue;

    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectJets(const cat::JetCollection& jets, const ParticleCollection& recolep) const
{
  cat::JetCollection seljets;
  for (auto& jet : jets) {
    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 2.4)  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep.p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}

cat::JetCollection TtbarDiLeptonAnalyzer::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.814) continue;
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
