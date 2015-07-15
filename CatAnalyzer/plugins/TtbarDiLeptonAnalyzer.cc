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
  vector<cat::Jet> selectJets(const edm::View<cat::Jet>* jets );
  vector<cat::Jet> selectBJets(vector<cat::Jet> & jets );

  TLorentzVector leafToTLorentzVector(reco::LeafCandidate & leaf)
  {return TLorentzVector(leaf.px(), leaf.py(),leaf.pz(),leaf.energy());}

  edm::EDGetTokenT<edm::View<cat::Muon> >     muonToken_;
  edm::EDGetTokenT<edm::View<cat::Electron> > elecToken_;
  edm::EDGetTokenT<edm::View<cat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<cat::MET> >      metToken_;
  edm::EDGetTokenT<reco::VertexCollection >   vtxToken_;

  TTree * ttree_;
  int b_nbjet;
  float b_MET, b_ll_mass, b_maxweight;

  TtFullLepKinSolver* solver;
  double tmassbegin_, tmassend_, tmassstep_;
  vector<double> nupars_;

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

  tmassbegin_     = iConfig.getParameter<double>       ("tmassbegin");
  tmassend_       = iConfig.getParameter<double>       ("tmassend");
  tmassstep_      = iConfig.getParameter<double>       ("tmassstep");
  nupars_         = iConfig.getParameter<vector<double> >("neutrino_parameters");
  
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("top", "top");
  ttree_->Branch("nbjet", &b_nbjet, "nbjet/I");
  ttree_->Branch("ll_mass", &b_ll_mass, "ll_mass/F");
  ttree_->Branch("MET", &b_MET, "MET/F");

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


  vector<cat::Muon> selectedMuons = selectMuons( muons.product() );
  vector<cat::Electron> selectedElectrons = selectElecs( electrons.product() );
  vector<cat::Jet> selectedJets = selectJets( jets.product() );
  vector<cat::Jet> selectedBJets = selectBJets( selectedJets );

  //  printf("selectedMuons %lu, selectedElectrons %lu, selectedJets %lu, selectedBJets %lu\n",selectedMuons.size(), selectedElectrons.size(), selectedJets.size(), selectedBJets.size() );

  vector<TLorentzVector> recolep; 
  for (auto lep : selectedMuons){ recolep.push_back(lep.tlv());}
  for (auto lep : selectedElectrons){recolep.push_back(lep.tlv());}

  if (recolep.size() < 2) return;
  
  TLorentzVector met = mets->front().tlv();

  b_ll_mass = 90.;  
  b_MET = met.Pt();    
  b_nbjet = selectedBJets.size();
  
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
  b_maxweight = maxweight;
  //  printf("maxweight %f, top1.M() %f, top2.M() %f \n",maxweight, top1.M(), top2.M() );

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
    //    printf("muon with pt %4.1f, POG loose id %d, tight id %d\n",mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }

  return selmuons;
}

vector<cat::Electron> TtbarDiLeptonAnalyzer::selectElecs(const edm::View<cat::Electron>* elecs )
{
  vector<cat::Electron> selelecs;
  for (auto el : *elecs) {
    if (el.pt() < 5) continue;
    selelecs.push_back(el);
  }
  return selelecs;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectJets(const edm::View<cat::Jet>* jets )
{
  vector<cat::Jet> seljets;
  for (auto jet : *jets) {
    if (jet.pt() < 5) continue;
    seljets.push_back(jet);
  }
  return seljets;
}

vector<cat::Jet> TtbarDiLeptonAnalyzer::selectBJets(vector<cat::Jet> & jets )
{
  vector<cat::Jet> selBjets;
  for (auto jet : jets) {
    if (jet.pt() < 5) continue;
    selBjets.push_back(jet);
  }
  return selBjets;
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
