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

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  h_nevents = fs->make<TH1D>("nevents","nevents",1,0,1);
  ttree_ = fs->make<TTree>("bjet_color","bjet_color");
  data = new Data(ttree_);

}

void JetChargeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);

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

  auto selectedJets  = selectJets(*jets, recolep );
  auto selectedBJets =  selectBJets( selectedJets );

  const bool runOnMC = !iEvent.isRealData();
  if ( !runOnMC  ) { std::cout<<"It is not MC samples"<<std::endl; exit(-1); }
  cat::JetCollection fromBJet; 
  for( auto& j : selectedJets ) {
    if ( abs( j.partonPdgId()) != 5) continue;
    fromBJet.push_back(j); 
  }
  while ( fromBJet.size()<2 ) {
    auto temp_jet = cat::Jet();
    temp_jet.setPartonPdgId(0);

    fromBJet.push_back( temp_jet );
    
  } 
  data->reset();
  for( int i=0 ; i<2 ; ++i) {
    data->lep_pt[i]= recolep[i]->pt();
    data->lep_pdgId[i] = recolep[i]->pdgId();
    data->jet_pt[i] = fromBJet[i].pt();
    data->jet_pdgId[i] = fromBJet[i].partonPdgId();
    data->jet_charge[i] = fromBJet[i].charge();
    LogDebug("JetChargeAnalyzer")<<"Jet Charge from b quark jet : "<<fromBJet[i].charge();
  }
  ttree_->Fill();
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
