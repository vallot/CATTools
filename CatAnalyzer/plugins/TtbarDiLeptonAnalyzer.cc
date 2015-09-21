#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"

#include "CATTools/CatAnalyzer/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

class TtbarDiLeptonAnalyzer : public edm::EDAnalyzer {
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
  vector<cat::Muon> selectMuons(const edm::View<cat::Muon>* muons );
  vector<cat::Electron> selectElecs(const edm::View<cat::Electron>* elecs );
  vector<cat::Jet> selectJets(const edm::View<cat::Jet>* jets, vector<cat::Particle> recolep);
  vector<cat::Jet> selectBJets(vector<cat::Jet> & jets );
  const reco::Candidate* getLast(const reco::Candidate* p);

  edm::EDGetTokenT<bool>          goodVertices_;
  edm::EDGetTokenT<bool>          CSCTightHaloFilter_;
  edm::EDGetTokenT<bool>          HBHENoiseFilter_;
  edm::EDGetTokenT<bool>          eeBadScFilter_;
  edm::EDGetTokenT<vector<pair<string, int> > > triggers_;

  edm::EDGetTokenT<edm::View<cat::Muon> >     muonToken_;
  edm::EDGetTokenT<edm::View<cat::Electron> > elecToken_;
  edm::EDGetTokenT<edm::View<cat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> >      metToken_;
  edm::EDGetTokenT<reco::VertexCollection >   vtxToken_;
  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<vector<int> > partonTop_modes_;
  edm::EDGetTokenT<reco::GenParticleCollection > partonTop_genParticles_;

  edm::EDGetTokenT<vector<reco::GenJet>      > pseudoTop_jets_;
  edm::EDGetTokenT<vector<reco::GenJet>      > pseudoTop_leptons_;
  edm::EDGetTokenT<vector<reco::GenParticle> > pseudoTop_;
  edm::EDGetTokenT<vector<reco::GenParticle> > pseudoTop_neutrinos_;
  edm::EDGetTokenT<vector<reco::MET>         > pseudoTop_mets_;

  TTree * ttree_;
  int b_partonChannel, b_partonMode1, b_partonMode2;
  int b_pseudoTopChannel, b_pseudoTopMode1, b_pseudoTopMode2;
  int b_njet, b_nbjet, b_step, b_channel, b_lepinPhase, b_jetinPhase;
  float b_MET, b_maxweight;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;
  float b_jet_pt, b_jet_eta, b_jet_phi, b_jet_m, b_jet_CSVInclV2;
  float b_top1_pt, b_top1_eta, b_top1_phi;
  float b_top2_pt, b_top2_eta, b_top2_phi;
  float b_tri;
  int b_filtered;
  int b_is3lep;

  std::unique_ptr<TtFullLepKinSolver> solver;
  bool runOnMC_;
  std::vector<int> pseudoTop_modes;
  //enum TTbarMode { CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON };
  //enum DecayMode { CH_HADRON = 1, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };
  int cutflow[10][4];
};
//
// constructors and destructor
//
TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
  goodVertices_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("goodVertices"));
  CSCTightHaloFilter_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("CSCTightHaloFilter"));
  HBHENoiseFilter_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("HBHENoiseFilter"));
  eeBadScFilter_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("eeBadScFilter"));
  triggers_      = consumes<vector<pair<string, int> > >(iConfig.getParameter<edm::InputTag>("triggers"));
  muonToken_ = consumes<edm::View<cat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<edm::View<cat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<edm::View<cat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<edm::View<cat::MET> >(iConfig.getParameter<edm::InputTag>("mets"));     
  vtxToken_  = consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  partonTop_channel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("partonTop_channel"));
  partonTop_modes_   = consumes<vector<int> >(iConfig.getParameter<edm::InputTag>("partonTop_modes"));
  partonTop_genParticles_   = consumes<reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("partonTop_genParticles"));
  pseudoTop_jets_      = consumes<vector<reco::GenJet>      >(iConfig.getParameter<edm::InputTag>("pseudoTop_jets"));
  pseudoTop_leptons_   = consumes<vector<reco::GenJet>      >(iConfig.getParameter<edm::InputTag>("pseudoTop_leptons"));
  pseudoTop_           = consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pseudoTop"));
  pseudoTop_neutrinos_ = consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pseudoTop_neutrinos"));
  pseudoTop_mets_      = consumes<vector<reco::MET>         >(iConfig.getParameter<edm::InputTag>("pseudoTop_mets"));

  const double tmassbegin = iConfig.getParameter<double>       ("tmassbegin");
  const double tmassend   = iConfig.getParameter<double>       ("tmassend");
  const double tmassstep  = iConfig.getParameter<double>       ("tmassstep");
  const auto   nupars     = iConfig.getParameter<vector<double> >("neutrino_parameters");
  
  solver.reset(new TtFullLepKinSolver(tmassbegin, tmassend, tmassstep, nupars));

  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("tree", "tree");
  ttree_->Branch("parton_channel", &b_partonChannel, "parton_channel/I");
  ttree_->Branch("parton_mode1", &b_partonMode1, "parton_mode1/I");
  ttree_->Branch("parton_mode2", &b_partonMode2, "parton_mode2/I");
  ttree_->Branch("pseudoTop_channel", &b_pseudoTopChannel, "pseudoTop_channel/I");
  ttree_->Branch("pseudoTop_mode1", &b_pseudoTopMode1, "pseudoTop_mode1/I");
  ttree_->Branch("pseudoTop_mode2", &b_pseudoTopMode2, "pseudoTop_mode2/I");

  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("nbjet", &b_nbjet, "nbjet/I");
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("step", &b_step, "step/I");
  ttree_->Branch("lepinPhase", &b_lepinPhase, "lepinPhase/I");
  ttree_->Branch("jetinPhase", &b_jetinPhase, "jetinPhase/I");

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
  ttree_->Branch("jet_pt", &b_jet_pt, "jet_pt/F");
  ttree_->Branch("jet_eta", &b_jet_eta, "jet_eta/F");
  ttree_->Branch("jet_phi", &b_jet_phi, "jet_phi/F");
  ttree_->Branch("jet_m", &b_jet_m, "jet_m/F");
  ttree_->Branch("jet_CSVInclV2", &b_jet_CSVInclV2, "jet_CSVInclV2/F");

  ttree_->Branch("top1_pt", &b_top1_pt, "top1_pt/F");
  ttree_->Branch("top1_eta", &b_top1_eta, "top1_eta/F");
  ttree_->Branch("top1_phi", &b_top1_pt, "top1_phi/F");
  ttree_->Branch("top2_pt", &b_top2_pt, "top2_pt/F");
  ttree_->Branch("top2_eta", &b_top2_eta, "top2_eta/F");
  ttree_->Branch("top2_phi", &b_top2_pt, "top2_phi/F");

  ttree_->Branch("tri", &b_tri, "tri/F");
  ttree_->Branch("filtered", &b_filtered, "filtered/I");
  ttree_->Branch("is3lep", &b_is3lep, "is3lep/I");

  for (int i = 0; i < 10; i++) for (int j = 0; j < 4; j++) cutflow[i][j]=0;
}

TtbarDiLeptonAnalyzer::~TtbarDiLeptonAnalyzer()
{
  cout <<"cut flow         emu         ee         mumu"<< endl;
  for (int i = 0; i < 6; i++){
    cout <<"step "<< i<< " "<< cutflow[i][0] <<  " "<< cutflow[i][1] << " " << cutflow[i][2] << " " << cutflow[i][3]<< endl;
  }
}

void TtbarDiLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_partonChannel = -1; b_partonMode1 = -1; b_partonMode2 = -1; 
  b_pseudoTopChannel = -1; b_pseudoTopMode1 = -1; b_pseudoTopMode2 = -1; 
  b_MET = -1; 
  b_njet = -1;
  b_nbjet = -1;
  b_channel = 0;
  b_step = -1;
  b_lepinPhase = 0; b_jetinPhase = 0;
  b_lep1_pt = -9; b_lep1_eta = -9; b_lep1_phi = -9;
  b_lep2_pt = -9; b_lep2_eta = -9; b_lep2_phi = -9;
  b_ll_pt = -9; b_ll_eta = -9; b_ll_phi = -9; b_ll_m = -9;
  b_jet_pt = -9; b_jet_eta = -9; b_jet_phi = -9; b_jet_m = -9; b_jet_CSVInclV2 = -9;
  b_top1_pt = -9; b_top1_eta = -9; b_top1_phi = -9;
  b_top2_pt = -9; b_top2_eta = -9; b_top2_phi = -9;
  b_tri = -9;
  b_filtered = -9; b_is3lep = -9;
  runOnMC_ = !iEvent.isRealData();

  // bool debug = false;
  // if (iEvent.id().event() == 312909020 || iEvent.id().event() == 255013550){
  //   debug = true;
  //   cout <<"############## debugging iEvent.id().event()" << iEvent.id().event()<<endl;
  // }
  edm::Handle<reco::VertexCollection> vertices;      iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ return;} // skip the event if no PV found
  // const reco::Vertex &PV = vertices->front();
  edm::Handle<edm::View<cat::Muon> > muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<edm::View<cat::Electron> > electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<edm::View<cat::Jet> > jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<edm::View<cat::MET> > mets;            iEvent.getByToken(metToken_, mets);
    
  if (runOnMC_){
    edm::Handle<int> partonTop_channel;
    edm::Handle<vector<int> > partonTop_modes;
    edm::Handle<reco::GenParticleCollection > partonTop_genParticles;    
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

    edm::Handle<vector<reco::GenJet>      > pseudoTop_jets;
    edm::Handle<vector<reco::GenJet>      > pseudoTop_leptons;
    edm::Handle<vector<reco::GenParticle> > pseudoTop;
    edm::Handle<vector<reco::GenParticle> > pseudoTop_neutrinos;
    edm::Handle<vector<reco::MET>         > pseudoTop_mets;
    iEvent.getByToken(pseudoTop_jets_     , pseudoTop_jets);
    iEvent.getByToken(pseudoTop_leptons_  , pseudoTop_leptons);
    iEvent.getByToken(pseudoTop_          , pseudoTop);
    iEvent.getByToken(pseudoTop_neutrinos_, pseudoTop_neutrinos);
    iEvent.getByToken(pseudoTop_mets_     , pseudoTop_mets);

    if ((*pseudoTop_leptons).size() == 2){
      if (((*pseudoTop_leptons)[0].pt() > 20) && ((*pseudoTop_leptons)[1].pt() > 20) && (std::abs((*pseudoTop_leptons)[0].eta()) < 2.4) && (std::abs((*pseudoTop_leptons)[1].eta()) < 2.4)) b_lepinPhase = 1;
    }

    if ((*pseudoTop_jets).size() == 2){
      if (((*pseudoTop_jets)[0].pt() > 30) && ((*pseudoTop_jets)[1].pt() > 30) && (std::abs((*pseudoTop_jets)[0].eta()) < 2.4) && (std::abs((*pseudoTop_jets)[1].eta()) < 2.4)) b_jetinPhase = 1;
    }

    int mode = 0;
    pseudoTop_modes.clear();
    for (const reco::GenJet & g : *pseudoTop_leptons){
      if ( std::abs(g.pdgId()) == 13){ mode = 2; }
      else if ( std::abs(g.pdgId()) == 11){ mode = 3; }
      pseudoTop_modes.push_back(mode);
    }

    if ( pseudoTop_modes.size() < 2 ){
      for (const reco::GenParticle & g : *pseudoTop_neutrinos){
	if (std::abs(g.pdgId()) == 16){
	  pseudoTop_modes.push_back(4);
	}
      }
    }

    if ( pseudoTop_modes.size() == 0 ) { pseudoTop_modes.push_back(0); }
    if ( pseudoTop_modes.size() == 1 ) { pseudoTop_modes.push_back(0); }
    b_pseudoTopMode1 = pseudoTop_modes[0]; 
    b_pseudoTopMode2 = pseudoTop_modes[1]; 

  }

  edm::Handle<bool> goodVertices;
  edm::Handle<bool> CSCTightHaloFilter;
  edm::Handle<bool> HBHENoiseFilter;
  edm::Handle<bool> eeBadScFilter;
  iEvent.getByToken(goodVertices_, goodVertices);
  iEvent.getByToken(CSCTightHaloFilter_, CSCTightHaloFilter);
  iEvent.getByToken(HBHENoiseFilter_, HBHENoiseFilter);
  iEvent.getByToken(eeBadScFilter_, eeBadScFilter);

  edm::Handle<vector<pair<string, int>>> triggers;
  iEvent.getByToken(triggers_, triggers);
  
  // if (!(*goodVertices && *CSCTightHaloFilter && *HBHENoiseFilter && *eeBadScFilter )){
  //   b_filtered = 1;
  //   return;
  // }
 
  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );
  vector<cat::Electron> selectedElectrons = selectElecs( electrons.product() );
  if (selectedMuons.size() + selectedElectrons.size() < 2) return;
    
  vector<cat::Particle> recolep; 
  for (auto lep : selectedMuons) recolep.push_back(lep);
  for (auto lep : selectedElectrons) recolep.push_back(lep);
  sort(recolep.begin(), recolep.end(), AnalysisHelper::ptSorting);
  
  int pdgIdSum = fabs(recolep[0].pdgId()) + fabs(recolep[1].pdgId());
  if (pdgIdSum == 24) b_channel = 1; // emu
  if (pdgIdSum == 22) b_channel = 2; // ee
  if (pdgIdSum == 26) b_channel = 3; // mumu

  bool tri=false;
  for (auto &t: *triggers){
    if (t.second){
      if (t.first.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") == 0 ||
    	t.first.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") == 0 )
      if (b_channel == 2) tri = true;
    
    if (t.first.find("HLT_Mu17_Mu8_DZ_v") == 0 ||
    	t.first.find("HLT_Mu17_TkMu8_DZ_v") == 0 ||	
    	t.first.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") == 0 ||	
    	t.first.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") == 0 ||
    	t.first.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") == 0 )
      if (b_channel == 3) tri = true;
    
    if (t.first.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") == 0 ||
	t.first.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v") == 0 )
      if (b_channel == 1) tri = true;
    //if (tri) cout << t.first << endl;
    }
  }
  if (!tri) return;
  
  cutflow[0][b_channel]++;

  b_lep1_pt = recolep[0].pt(); b_lep1_eta = recolep[0].eta(); b_lep1_phi = recolep[0].phi();
  b_lep2_pt = recolep[1].pt(); b_lep2_eta = recolep[1].eta(); b_lep2_phi = recolep[1].phi();
  TLorentzVector tlv_ll = recolep[0].tlv()+recolep[1].tlv();
  b_ll_pt = tlv_ll.Pt(); b_ll_eta = tlv_ll.Eta(); b_ll_phi = tlv_ll.Phi(); b_ll_m = tlv_ll.M();

  if (b_ll_m < 20.){
    ttree_->Fill();
    return;
  }
  if (recolep[0].charge() * recolep[1].charge() > 0){
    ttree_->Fill();
    return;
  }
  cutflow[1][b_channel]++;
  b_step = 1;
 
  if (b_channel != 1){
    if ((b_ll_m > 76) && (b_ll_m < 106)){
      ttree_->Fill();
      return;
    }
  }
  cutflow[2][b_channel]++;
  b_step = 2;

  vector<cat::Jet> selectedJets = selectJets( jets.product(), recolep );
  vector<cat::Jet> selectedBJets = selectBJets( selectedJets );
  TLorentzVector met = mets->front().tlv();
  b_MET = met.Pt();
  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();

  if (selectedJets.size() < 2){
    ttree_->Fill();
    return;
  }
  cutflow[3][b_channel]++;
  b_step = 3;

  if (b_channel != 1){
    if (b_MET < 40.){
      ttree_->Fill();
      return;
    }
  }
  cutflow[4][b_channel]++;
  b_step = 4;
  
  if (selectedBJets.size() == 0){
    ttree_->Fill();
    return;
  }
  cutflow[5][b_channel]++;
  b_step = 5;

  //  printf("selectedMuons %lu, selectedElectrons %lu, selectedJets %lu, selectedBJets %lu\n",selectedMuons.size(), selectedElectrons.size(), selectedJets.size(), selectedBJets.size() );
  /*
    for (auto jet : selectedJets) {
    TLorentzVector tlv_jet = jet.tlv();
    b_jet_pt = tlv_jet.Pt();
    b_jet_eta = tlv_jet.Eta();
    b_jet_phi = tlv_jet.Phi();
    b_jet_m = tlv_jet.M();
    b_jet_CSVInclV2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    }
  */

  ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
  int kin=0; TLorentzVector nu1, nu2, top1, top2;
  double maxweight=0;
  cat::Jet kinj1, kinj2;

  for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
    for (auto jet2 = next(jet1); jet2 != end; ++jet2){

      double weight1 =0; double weight2 =0;
      TLorentzVector recojet1= jet1->tlv();
      TLorentzVector recojet2= jet2->tlv();
      TLorentzVector recolep1= recolep[0].tlv();
      TLorentzVector recolep2= recolep[1].tlv();
      
      double xconstraint = recolep1.Px()+recolep2.Px()+ (recojet1).Px() + (recojet2).Px() +met.Px();
      double yconstraint = recolep1.Py()+recolep2.Py()+ (recojet2).Py() + (recojet1).Py() +met.Py();
      
      solver->SetConstraints(xconstraint, yconstraint);
      TtFullLepKinSolver::NeutrinoSolution nuSol= solver->getNuSolution( recolep1, recolep2 , recojet1, recojet2);
      weight1 = nuSol.weight;
      TLorentzVector nu11 = AnalysisHelper::leafToTLorentzVector(nuSol.neutrino);
      TLorentzVector nu12 = AnalysisHelper::leafToTLorentzVector(nuSol.neutrinoBar);
      
      TtFullLepKinSolver::NeutrinoSolution nuSol2= solver->getNuSolution( recolep1, recolep2 , recojet2, recojet1);
      weight2 = nuSol2.weight;
      TLorentzVector nu21 = AnalysisHelper::leafToTLorentzVector(nuSol2.neutrino);
      TLorentzVector nu22 = AnalysisHelper::leafToTLorentzVector(nuSol2.neutrinoBar);
      if (weight1 > maxweight || weight2 > maxweight){
        if(weight1>weight2 && weight1>0){
          maxweight = weight1; kinj1=(*jet1); kinj2=(*jet2); nu1 = nu11; nu2 = nu12; kin++;
          top1 = recolep1+recojet1+nu11; top2 = recolep2+recojet2+nu12;
        }
        else if(weight2>weight1 && weight2>0){
          maxweight = weight2; kinj1=(*jet2); kinj2=(*jet1); nu1 = nu21; nu2 = nu22; kin++;
          top1 = recolep1+recojet2+nu21; top2 = recolep2+recojet1+nu22;
        }
      }
    }
  }

  b_top1_pt = top1.Pt();
  b_top1_eta = top1.Eta();
  b_top1_phi = top1.Phi();
  b_top2_pt = top2.Pt();
  b_top2_eta = top2.Eta();
  b_top2_phi = top2.Phi();

  b_maxweight = maxweight;
  //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
  // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_MET, b_ll_mass, b_maxweight);

  ttree_->Fill();
}

const reco::Candidate* TtbarDiLeptonAnalyzer::getLast(const reco::Candidate* p)
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
    {
      const reco::Candidate* dau = p->daughter(i);
      if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
    }
  return p;
}

vector<cat::Muon> TtbarDiLeptonAnalyzer::selectMuons(const edm::View<cat::Muon>* muons )
{
  vector<cat::Muon> selmuons;
  for (auto mu : *muons) {
    if (mu.pt() < 20.) continue;
    if (fabs(mu.eta()) > 2.4) continue;
    //if (!mu.isMediumMuon()) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.12) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }

  return selmuons;
}

vector<cat::Electron> TtbarDiLeptonAnalyzer::selectElecs(const edm::View<cat::Electron>* elecs )
{
  vector<cat::Electron> selelecs;
  for (auto el : *elecs) {
    if (el.pt() < 20.) continue;
    if ((fabs(el.scEta()) > 1.4442) && (fabs(el.scEta()) < 1.566)) continue;
    if (fabs(el.eta()) > 2.4) continue;
    //if (!el.electronID("cutBasedElectronID-Spring15-50ns-V1-standalone-medium")) continue;
    if (el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium") == 0) continue;
    //if (!el.passConversionVeto()) continue;
    if (!el.isPF()) continue;
        
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return selelecs;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectJets(const edm::View<cat::Jet>* jets, vector<cat::Particle> recolep )
{
  vector<cat::Jet> seljets;seljets.clear();
  for (auto jet : *jets) {
    if (jet.pt() < 30.) continue;
    if (fabs(jet.eta()) > 2.4)	continue;
    if (!jet.LooseId()) continue;
    
    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet,lep) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectBJets(vector<cat::Jet> & jets )
{
  vector<cat::Jet> selBjets;
  for (auto jet : jets) {
    if (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.605) continue;	
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TtbarDiLeptonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
