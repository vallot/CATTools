#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/Jet.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CatAnalyzer/plugins/JetChargeAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TTree.h"
#include "TH1D.h"

using namespace std;
using namespace cat;

//
// constructors and destructor
//




JetChargeAnalyzer::JetChargeAnalyzer(const edm::ParameterSet& iConfig)
{

  const auto muonSet = iConfig.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  const auto elecSet = iConfig.getParameter<edm::ParameterSet>("electron");
  elecToken_ = consumes<cat::ElectronCollection>(elecSet.getParameter<edm::InputTag>("src"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  solverCandsToken_  = consumes<std::vector<reco::LeafCandidate> >(iConfig.getParameter<edm::InputTag>("solverCands"));
  solverQualityToken_  = consumes<std::vector<float> >(iConfig.getParameter<edm::InputTag>("solverQuality"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  h_nevents = fs->make<TH1D>("nevents","nevents",1,0,1);
  ttree_ = fs->make<TTree>("true_bjet_charge","true_bjet_charge");
  rtree_ = fs->make<TTree>("reco_bjet_charge","reco_bjet_charge");
  data1 = new Data(ttree_);
  data2 = new Data(rtree_);

  //auto solverPSet = iConfig.getParameter<edm::ParameterSet>("solver");
  //solver_.reset( new CMSKinSolver(solverPSET));  


}

void JetChargeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  edm::Handle<cat::MuonCollection> muons;                         iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;                 iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;                           iEvent.getByToken(jetToken_, jets);
  edm::Handle<std::vector<reco::LeafCandidate> > solverCands;     iEvent.getByToken(solverCandsToken_, solverCands);
  edm::Handle<std::vector<float> > qualities;                    iEvent.getByToken(solverQualityToken_, qualities);
  // Find leptons and sort by pT
  cat::MuonCollection selMuons;
  cat::ElectronCollection selElecs;

  selectMuons(*muons, selMuons);
  selectElecs(*electrons, selElecs);
  if ( selMuons.size()+selElecs.size() < 2 ) {
    return ;
  }

  std::vector<const cat::Lepton*> recolep;
  for ( const auto& x : selMuons ) recolep.push_back(&x);
  for ( const auto& x : selElecs ) recolep.push_back(&x);

  sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
  h_nevents->Fill(1);
  auto selectedJets  = selectJets(*jets, recolep );
  auto selectedBJets =  selectBJets( selectedJets );

  const bool runOnMC = !iEvent.isRealData();
  if ( !runOnMC  ) { std::cout<<"It is not MC samples"<<std::endl; exit(-1); }
  cat::JetCollection fromTrueBJet, fromRecoBJet; 

  // First b-related jets
  for( auto& j : selectedJets ) {
    if ( abs( j.partonPdgId()) != 5) continue;
    fromTrueBJet.push_back(j); 
  }
  for( auto& j : selectedJets ) {
    if ( !j.CSVv2L() ) continue;
    fromRecoBJet.push_back(j); 
  }
  // Second : Normal jet( no b-releated jets )
  if ( fromTrueBJet.size()<2 ) {
    for( auto& j : selectedJets) {
      if ( abs(j.partonPdgId()) ==5 ) continue;
      fromTrueBJet.push_back( j );
    }
  } 
  if ( fromRecoBJet.size()<2 ) {
    for( auto& j : selectedJets) {
      if ( j.CSVv2L() ) continue;
      fromRecoBJet.push_back( j );
    }
  } 
  // Last : Dummy jet
  while ( fromTrueBJet.size()<2 ) {
    auto temp_jet = cat::Jet();
    temp_jet.setPartonPdgId(0);
    fromTrueBJet.push_back( temp_jet );
  } 
  while ( fromRecoBJet.size()<2 ) {
    auto temp_jet = cat::Jet();
    temp_jet.setPartonPdgId(0);
    fromRecoBJet.push_back( temp_jet );
  }
  // However, only first 2 jets wll be used. 
  data1->reset();
  for( int i=0 ; i<2 ; ++i) {
    data1->lep_pt[i]= recolep[i]->pt();
    data1->lep_pdgId[i] = recolep[i]->pdgId();
    data1->jet_pt[i] = fromTrueBJet[i].pt();
    data1->jet_pdgId[i] = fromTrueBJet[i].partonPdgId();
    if ( data1->jet_pdgId[i] ==0 ) {
      data1->jet_eta[i] = -999.f;
      data1->jet_charge[i] = -999;
    }
    else {
      data1->jet_eta[i] = fromTrueBJet[i].eta();
      data1->jet_charge[i] = fromTrueBJet[i].charge();
    }
    data1->jet_btag[i] = fromTrueBJet[i].CSVv2L();
    LogDebug("JetChargeAnalyzer")<<"Jet Charge from true b quark jet : "<<fromTrueBJet[i].charge();
  }
  data1->solverQuality = (float)((*qualities)[0]*1e5);
  data1->top_mass[0] = (*solverCands)[1].mass();
  data1->top_mass[1] = (*solverCands)[2].mass();
  auto bjet1 = (*solverCands)[7].p4();
  auto bjet2 = (*solverCands)[8].p4();

  data1->bjet[0] = new TLorentzVector(bjet1.px(), bjet1.py(), bjet1.pz(), bjet1.energy());
  data1->bjet[1] = new TLorentzVector(bjet2.px(), bjet2.py(), bjet2.pz(), bjet2.energy());
  data1->bjet_charge[0] = (*solverCands)[7].charge();
  data1->bjet_charge[1] = (*solverCands)[8].charge();

  data2->reset();
  for( int i=0 ; i<2 ; ++i) {
    data2->lep_pt[i]= recolep[i]->pt();
    data2->lep_pdgId[i] = recolep[i]->pdgId();
    data2->jet_pt[i] = fromRecoBJet[i].pt();
    //data2->jet_eta[i] = fromRecoBJet[i].eta();
    data2->jet_pdgId[i] = fromRecoBJet[i].partonPdgId();
    if ( data2->jet_pdgId[i] ==0 ) {
      data2->jet_eta[i] = -999.f;
      data2->jet_charge[i] = -999;
    }
    else {
      data2->jet_eta[i] = fromRecoBJet[i].eta();
      data2->jet_charge[i] = fromRecoBJet[i].charge();
    }
    data2->jet_btag[i] = fromRecoBJet[i].CSVv2L();
    LogDebug("JetChargeAnalyzer")<<"Jet Charge from reco b quark jet : "<<fromRecoBJet[i].charge();
  }
  data2->solverQuality = (float)((*qualities)[0]*1e5);
  data2->top_mass[0] = (*solverCands)[1].mass();
  data2->top_mass[1] = (*solverCands)[2].mass();

  data2->bjet[0] = new TLorentzVector(bjet1.px(), bjet1.py(), bjet1.pz(), bjet1.energy());
  data2->bjet[1] = new TLorentzVector(bjet2.px(), bjet2.py(), bjet2.pz(), bjet2.energy());
  data2->bjet_charge[0] = (*solverCands)[7].charge();
  data2->bjet_charge[1] = (*solverCands)[8].charge();

  ttree_->Fill();
  rtree_->Fill();
}


float JetChargeAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons ) const
{
  float weight = 1.;
  for (auto& m : muons) {
    cat::Muon mu(m);

    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.15) continue;
    selmuons.push_back(mu);
  }
  return weight;
}

float JetChargeAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs) const
{
  float weight = 1.;
  for (auto& e : elecs) {
    cat::Electron el(e);

    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) continue;
    if (el.relIso(0.3) > 0.12) continue;

    selelecs.push_back(el);
  }
  return weight;
}
cat::JetCollection JetChargeAnalyzer::selectJets(const cat::JetCollection& jets, const std::vector<const cat::Lepton*>& recolep)
{
  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);

    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 2.4)  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep->p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    seljets.push_back(jet);
  }
  return seljets;
}

cat::JetCollection JetChargeAnalyzer::selectBJets(const JetCollection& jets) const
{
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2L) continue;
    //if (jet.bDiscriminator(BTAG_CSVv2) < WP_BTAG_CSVv2M) continue;//forsync
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetChargeAnalyzer);
