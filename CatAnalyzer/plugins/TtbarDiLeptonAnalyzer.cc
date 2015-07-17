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

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

class TtbarDiLeptonAnalyzer : public edm::EDAnalyzer {
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  vector<cat::Muon> selectMuons(const edm::View<cat::Muon>* muons );
  vector<cat::Electron> selectElecs(const edm::View<cat::Electron>* elecs );
  vector<cat::Jet> selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep);
  vector<cat::Jet> selectBJets(vector<cat::Jet> & jets );
  float passingSteps(int channel, float met, float ll_mass, float ll_charge, int selectedJets_size, int btag);

  TLorentzVector leafToTLorentzVector(reco::LeafCandidate & leaf)
  {return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}

  edm::EDGetTokenT<edm::View<cat::Muon> >     muonToken_;
  edm::EDGetTokenT<edm::View<cat::Electron> > elecToken_;
  edm::EDGetTokenT<edm::View<cat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> >      metToken_;
  edm::EDGetTokenT<reco::VertexCollection >   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<vector<int> > partonTop_modes_;

  TTree * ttree_;
  int b_njet, b_nbjet, b_step, b_channel;
  float b_MET, b_ll_mass, b_maxweight;

  TtFullLepKinSolver* solver;
  double tmassbegin_, tmassend_, tmassstep_;
  vector<double> nupars_;
  
  bool runOnMC_;
};
//
// constructors and destructor
//
TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
  muonToken_ = consumes<edm::View<cat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<edm::View<cat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<edm::View<cat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<edm::View<cat::MET> >(iConfig.getParameter<edm::InputTag>("mets"));     
  vtxToken_  = consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));
  partonTop_channel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("partonTop_channel"));
  partonTop_modes_   = consumes<vector<int> >(iConfig.getParameter<edm::InputTag>("partonTop_modes"));

  tmassbegin_     = iConfig.getParameter<double>       ("tmassbegin");
  tmassend_       = iConfig.getParameter<double>       ("tmassend");
  tmassstep_      = iConfig.getParameter<double>       ("tmassstep");
  nupars_         = iConfig.getParameter<vector<double> >("neutrino_parameters");
  
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("top", "top");
  ttree_->Branch("njet", &b_njet, "njet/I");
  ttree_->Branch("nbjet", &b_nbjet, "nbjet/I");
  ttree_->Branch("ll_mass", &b_ll_mass, "ll_mass/F");
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("step", &b_step, "njet/I");

}
TtbarDiLeptonAnalyzer::~TtbarDiLeptonAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void TtbarDiLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle<edm::View<cat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<edm::View<cat::Electron> > electrons;
  iEvent.getByToken(elecToken_, electrons);

  edm::Handle<edm::View<cat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<edm::View<cat::MET> > mets;
  iEvent.getByToken(metToken_, mets);
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<int> partonTop_channel;
  edm::Handle<vector<int> > partonTop_modes;

  if (runOnMC_){
    iEvent.getByToken(mcLabel_,genParticles);
    /*
    for (auto g : *genParticles)
    {
      if (g.pdgId() == 6)
      {
        if (g.numberOfDaughters() == 0) continue;
//        const reco::GenParticle* partonTop1 = &g;
//        const reco::Candidate* p = partonTop1->daughter(0);
        const reco::Candidate* p = &g->daughter(0);
//        cout <<"partonW1->pdgId() " << partonW1->pdgId() <<endl; 
//        cout <<"g.pdgId() " << g.pdgId() <<endl; 
//        cout <<"there's g.daughter" <<endl; 
//        auto p = g.daughter(0);
//        reco::Candidate* p = g.daughter(0);
 //       cout <<"p.pdgId() " << p->pdgId() <<endl; 
        bool isLast = true;
        while (isLast)
        {
          if (p->pdgId() != p->daughter(0)->pdgId()) isLast = false;
          p = p->daughter(0);
        }
        cout <<"p.pdgId() " << p->pdgId() <<endl; 
        cout <<"p.daughter().pdgId() " << p->daughter(0)->pdgId() <<endl; 
      }
    }
    */
    iEvent.getByToken(partonTop_channel_, partonTop_channel);
    iEvent.getByToken(partonTop_modes_, partonTop_modes);
    cout << "partonTop_channel "<< *partonTop_channel<<endl;
    for (auto m : *partonTop_modes){
      cout << "partonTop_mode "<< m<<endl;
    }
	
  }
  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );
  vector<cat::Electron> selectedElectrons = selectElecs( electrons.product() );

  vector<TLorentzVector> recolep; 
  for (auto lep : selectedMuons){ recolep.push_back(lep.tlv()); }
  for (auto lep : selectedElectrons){ recolep.push_back(lep.tlv()); }
  if (recolep.size() < 2) return;

  float channel = selectedElectrons.size();

  float ll_charge = 0. ;
  if (channel == 0) ll_charge = selectedMuons[0].charge()*selectedMuons[1].charge();
  if (channel == 1) ll_charge = selectedMuons[0].charge()*selectedElectrons[0].charge();
  if (channel == 2) ll_charge = selectedElectrons[0].charge()*selectedElectrons[1].charge();

  vector<cat::Jet> selectedJets = selectJets( jets.product(), recolep );
  vector<cat::Jet> selectedBJets = selectBJets( selectedJets );

  //  printf("selectedMuons %lu, selectedElectrons %lu, selectedJets %lu, selectedBJets %lu\n",selectedMuons.size(), selectedElectrons.size(), selectedJets.size(), selectedBJets.size() );

  TLorentzVector met = mets->front().tlv();

  b_ll_mass = (recolep[0]+recolep[1]).M();  
  b_MET = met.Pt(); 
  b_njet = selectedJets.size();
  b_nbjet = selectedBJets.size();
  b_channel = channel;

  float step = passingSteps( channel, met.Pt(), (recolep[0]+recolep[1]).M(), ll_charge, selectedJets.size(), selectedBJets.size() );
  b_step = step;
  
  ////////////////////////////////////////////////////////  KIN  /////////////////////////////////////
  int kin=0; TLorentzVector nu1, nu2, top1, top2;
  double maxweight=0;
  cat::Jet kinj1, kinj2;

  for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
    for (auto jet2 = next(jet1); jet2 != end; ++jet2){

      double weight1 =0; double weight2 =0;
      TLorentzVector nu11, nu12, nu21, nu22;
      TLorentzVector recojet1= jet1->tlv();
      TLorentzVector recojet2= jet2->tlv();
      
      double xconstraint = recolep[0].Px()+recolep[1].Px()+ (recojet1).Px() + (recojet2).Px() +met.Px();
      double yconstraint = recolep[0].Py()+recolep[1].Py()+ (recojet2).Py() + (recojet1).Py() +met.Py();
      
      solver->SetConstraints(xconstraint, yconstraint);
      TtFullLepKinSolver::NeutrinoSolution nuSol= solver->getNuSolution( recolep[0], recolep[1] , recojet1, recojet2);
      weight1 = nuSol.weight;
      nu11 = leafToTLorentzVector(nuSol.neutrino);
      nu12 = leafToTLorentzVector(nuSol.neutrinoBar);
      
      TtFullLepKinSolver::NeutrinoSolution nuSol2= solver->getNuSolution( recolep[0], recolep[1] , recojet2, recojet1);
      weight2 = nuSol2.weight;
      nu21 = leafToTLorentzVector(nuSol2.neutrino);
      nu22 = leafToTLorentzVector(nuSol2.neutrinoBar);
      if (weight1 > maxweight || weight2 > maxweight){
	if(weight1>weight2 && weight1>0){
	  maxweight = weight1; kinj1=(*jet1); kinj2=(*jet2); nu1 = nu11; nu2 = nu12; kin++;
	  top1 = recolep[0]+recojet1+nu11; top2 = recolep[1]+recojet2+nu12;
	}
	else if(weight2>weight1 && weight2>0){
	  maxweight = weight2; kinj1=(*jet2); kinj2=(*jet1); nu1 = nu21; nu2 = nu22; kin++;
	  top1 = recolep[0]+recojet2+nu21; top2 = recolep[1]+recojet1+nu22;
	}
      }
    }
  }
  b_maxweight = maxweight;
  //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );
 // printf("%2d, %2d, %2d, %2d, %6.2f, %6.2f, %6.2f\n", b_njet, b_nbjet, b_step, b_channel, b_MET, b_ll_mass, b_maxweight);

  ttree_->Fill();
}

vector<cat::Muon> TtbarDiLeptonAnalyzer::selectMuons(const edm::View<cat::Muon>* muons )
{
  vector<cat::Muon> selmuons;
  for (auto mu : *muons) {
    if (!mu.isTightMuon()) continue;
    if (mu.pt() <= 20.) continue;
    if (fabs(mu.eta()) >= 2.4) continue;
    if (mu.relIso(0.4) >= 0.12) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }

  return selmuons;
}

vector<cat::Electron> TtbarDiLeptonAnalyzer::selectElecs(const edm::View<cat::Electron>* elecs )
{
  vector<cat::Electron> selelecs;
  for (auto el : *elecs) {
    if (!el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium")) continue;
    if (!el.passConversionVeto()) continue;
    if (!el.isPF()) continue;
    if (el.pt() <= 20.) continue;
    if ((fabs(el.scEta()) <= 1.4442) && (el.relIso(0.3) >= 0.1649)) continue;
    if ((fabs(el.scEta()) >= 1.566) && (el.relIso(0.3) >= 0.2075)) continue;
    if ((fabs(el.scEta()) > 1.4442) && (fabs(el.scEta()) < 1.566)) continue;
    if (fabs(el.eta()) >= 2.5) continue;
    if (el.pt() < 5) continue;
    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
  return selelecs;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectJets(const edm::View<cat::Jet>* jets, vector<TLorentzVector> recolep )
{
  vector<cat::Jet> seljets;
  for (auto jet : *jets) {
    if (!jet.LooseId()) continue;
    if (jet.pt() <= 30.) continue;
    if (fabs(jet.eta()) >= 2.4)	continue;
    if (jet.tlv().DeltaR(recolep[0]) <= 0.4) continue;
    if (jet.tlv().DeltaR(recolep[1]) <= 0.4) continue;
    // printf("jet with pt %4.1f\n", jet.pt());
    seljets.push_back(jet);
  }
  return seljets;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectBJets(vector<cat::Jet> & jets )
{
  vector<cat::Jet> selBjets;
  for (auto jet : jets) {
    float jets_CSVInclV2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    if (jets_CSVInclV2 <= 0.814) continue;	
    //printf("b jet with pt %4.1f\n", jet.pt());
    selBjets.push_back(jet);
  }
  return selBjets;
}

float TtbarDiLeptonAnalyzer::passingSteps(int channel, float met, float ll_mass, float ll_charge, int selectedJets_size, int btag)
{
  int step = 0;
  if (ll_mass <= 20.) return step;
  if (ll_charge > 0.) return step;
  step = 1;
  if (channel != 1){
    if ((ll_mass > 76) and (ll_mass < 106)) return step;
  }
  step = 2;
  if (selectedJets_size < 2) return step;
  step = 3;
  if (channel == 1){
    step = 4;
  }
  else{
    if (met <= 40.) return step;
  }
  step = 4;
  if (btag <= 0) return step;
  step = 5;

  return step;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TtbarDiLeptonAnalyzer::beginJob()
{
  solver = new TtFullLepKinSolver(tmassbegin_, tmassend_, tmassstep_, nupars_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TtbarDiLeptonAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TtbarDiLeptonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TtbarDiLeptonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TtbarDiLeptonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TtbarDiLeptonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
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
