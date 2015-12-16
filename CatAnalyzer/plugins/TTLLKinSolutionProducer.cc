#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "CATTools/CatAnalyzer/interface/KinematicSolvers.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace std;

namespace cat {

class TTLLKinSolutionProducer : public edm::stream::EDProducer<>
{
public:
  TTLLKinSolutionProducer(const edm::ParameterSet& pset);
  virtual ~TTLLKinSolutionProducer() {};
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  void produce(edm::Event & event, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<edm::View<reco::CandidatePtr> > leptonPtrToken_;
  edm::EDGetTokenT<edm::View<reco::CandidatePtr> > jetPtrToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > leptonToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > jetToken_;
  edm::EDGetTokenT<float> metToken_, metphiToken_;

private:
  typedef reco::Candidate::LorentzVector LV;
  typedef reco::LeafCandidate Cand;
  typedef std::vector<Cand> CandColl;
  typedef std::vector<float> floats;
  std::unique_ptr<KinematicSolver> solver_;

};

}

using namespace cat;

TTLLKinSolutionProducer::TTLLKinSolutionProducer(const edm::ParameterSet& pset)
{
  leptonPtrToken_ = mayConsume<edm::View<reco::CandidatePtr> >(pset.getParameter<edm::InputTag>("leptons"));
  leptonToken_ = mayConsume<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("leptons"));
  jetPtrToken_ = mayConsume<edm::View<reco::CandidatePtr> >(pset.getParameter<edm::InputTag>("jets"));
  jetToken_ = mayConsume<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("jets"));
  metToken_ = consumes<float>(pset.getParameter<edm::InputTag>("met"));
  metphiToken_ = consumes<float>(pset.getParameter<edm::InputTag>("metphi"));

  auto solverPSet = pset.getParameter<edm::ParameterSet>("solver");
  auto algoName = solverPSet.getParameter<std::string>("algo");
  std::transform(algoName.begin(), algoName.end(), algoName.begin(), ::toupper);
  if      ( algoName == "CMSKIN" ) solver_.reset(new CMSKinSolver(solverPSet));
  else if ( algoName == "DESYMASSLOOP" ) solver_.reset(new DESYMassLoopSolver(solverPSet));
  else if ( algoName == "DESYSMEARED" ) solver_.reset(new DESYSmearedSolver(solverPSet));
  else if ( algoName == "MT2"    ) solver_.reset(new MT2Solver(solverPSet));
  else if ( algoName == "MAOS"   ) solver_.reset(new MAOSSolver(solverPSet));
  else if ( algoName == "NUWGT"  ) solver_.reset(new NuWeightSolver(solverPSet));
  else if ( algoName == "DEFAULT" ) solver_.reset(new TTDileptonSolver(solverPSet));
  else {
    cerr << "The solver name \"" << solverPSet.getParameter<std::string>("algo") << "\" is not known please check spellings.\n";
    cerr << "Fall back to the default dummy solver\n";
    solver_.reset(new TTDileptonSolver(solverPSet)); // A dummy solver
  }

  produces<CandColl>();

  produces<floats>("aux");
  produces<floats>("mLL");
  produces<floats>("mLB");
  produces<floats>("mAddJJ");
  produces<floats>("dphi");
}

void TTLLKinSolutionProducer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&)
{
  if ( dynamic_cast<DESYSmearedSolver*>(solver_.get()) != 0 ) {
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(lumi.index());
    dynamic_cast<DESYSmearedSolver*>(solver_.get())->setRandom(&engine);
  }
}

void TTLLKinSolutionProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  std::auto_ptr<CandColl> cands(new CandColl);
  //auto candsRefProd = event.getRefBeforePut<CRCandColl>();

  std::auto_ptr<floats> out_aux(new floats);
  std::auto_ptr<floats> out_mLL(new floats);
  std::auto_ptr<floats> out_mLB(new floats);
  std::auto_ptr<floats> out_mAddJJ(new floats);
  std::auto_ptr<floats> out_dphi(new floats);

  std::vector<reco::CandidatePtr> leptons;
  edm::Handle<edm::View<reco::CandidatePtr> > leptonPtrHandle;
  edm::Handle<edm::View<reco::Candidate> > leptonHandle;
  if  ( event.getByToken(leptonPtrToken_, leptonPtrHandle) ) {
    for ( auto x : *leptonPtrHandle ) leptons.push_back(x);
  }
  else {
    event.getByToken(leptonToken_, leptonHandle);
    for ( int i=0, n=leptonHandle->size(); i<n; ++i ) leptons.push_back(reco::CandidatePtr(leptonHandle, i));
  }

  std::vector<reco::CandidatePtr> jets;
  edm::Handle<edm::View<reco::CandidatePtr> > jetPtrHandle;
  edm::Handle<edm::View<reco::Candidate> > jetHandle;
  if ( event.getByToken(jetPtrToken_, jetPtrHandle) ) {
    for ( auto x : *jetPtrHandle ) jets.push_back(x);
  }
  else {
    event.getByToken(jetToken_, jetHandle);
    for ( int i=0, n=jetHandle->size(); i<n; ++i ) jets.push_back(reco::CandidatePtr(jetHandle, i));
  }

  edm::Handle<float> metHandle;
  event.getByToken(metToken_, metHandle);
  const float met = *metHandle;
  event.getByToken(metphiToken_, metHandle);
  const float metphi = *metHandle;
  const LV metLV(met*cos(metphi), met*sin(metphi), 0, met);

  do {
    // Check objects to exist
    if ( leptons.size() < 2 ) break;
    if ( jets.size() < 2 ) break;

    // Pick leading leptons.
    const auto lep1 = leptons.at(0);
    const auto lep2 = leptons.at(1);
    const LV lep1LV = lep1->p4();
    const LV lep2LV = lep2->p4();
    LV inputLV[5] = {metLV, lep1LV, lep2LV};
    LV nu1LV, nu2LV;
    double quality = -1e9; // Default quality value

    // Run the solver with all jet combinations
    reco::CandidatePtr selectedJet1, selectedJet2;
    for ( auto jet1 : jets )
    {
      inputLV[3] = jet1->p4();
      for ( auto jet2 : jets )
      {
        if ( jet1 == jet2 ) continue;
        inputLV[4] = jet2->p4();

        solver_->solve(inputLV);
        if ( solver_->quality() > quality )
        {
          quality = solver_->quality();
          selectedJet1 = jet1;
          selectedJet2 = jet2;
          nu1LV = solver_->nu1();
          nu2LV = solver_->nu2();
        }
      }
    }
    if ( quality <= -1e9 ) break; // failed to get solution

    // Redo the calculation with the selected ones to update internal variables
    inputLV[3] = selectedJet1->p4();
    inputLV[4] = selectedJet2->p4();
    solver_->solve(inputLV);
    quality = solver_->quality();
    if ( quality <= -1e9 ) break;
    nu1LV = solver_->nu1();
    nu2LV = solver_->nu2();
    std::copy(solver_->aux().begin(), solver_->aux().end(), std::back_inserter(*out_aux));

    cands->resize(7);

    Cand& ttbar = cands->at(0);
    Cand& top1 = cands->at(1);
    Cand& top2 = cands->at(2);
    Cand& w1 = cands->at(3);
    Cand& w2 = cands->at(4);
    Cand& nu1 = cands->at(5);
    Cand& nu2 = cands->at(6);

    // Set four momentum
    nu1.setP4(nu1LV);
    nu2.setP4(nu2LV);
    w1.setP4(solver_->l1()+nu1LV);
    w2.setP4(solver_->l2()+nu2LV);
    top1.setP4(w1.p4()+solver_->j1());
    top2.setP4(w2.p4()+solver_->j2());
    ttbar.setP4(top1.p4()+top2.p4());

    // Set basic quantum numbers (do channel dependent things later)
    const int lep1Q = lep1->charge();
    const int lep2Q = lep2->charge();
    ttbar.setPdgId(0);
    w1.setPdgId(-24*lep1Q);
    w2.setPdgId(-24*lep2Q);
    top1.setPdgId(-6*lep1Q);
    top2.setPdgId(-6*lep2Q);

    // Do the basic mother-daughter associations
    /*
    auto prodId = candsRefProd.id();
    auto getter = candsRefProd.productGetter();
    ttbar.addDaughter(reco::CandidatePtr(prodId, 1, getter));
    ttbar.addDaughter(reco::CandidatePtr(prodId, 2, getter));
    top1.addDaughter(reco::CandidatePtr(prodId, 3, getter));
    top2.addDaughter(reco::CandidatePtr(prodId, 4, getter));

    nu1.setPdgId((pdgId1+1)*lep1->charge());
    nu2.setPdgId((pdgId2+1)*lep2->charge());

    w1.addDaughter(reco::CandidatePtr(prodId, 5, getter));
    w2.addDaughter(reco::CandidatePtr(prodId, 6, getter));
    */

    out_mLL->push_back((solver_->l1()+solver_->l2()).mass());
    out_dphi->push_back(deltaPhi(top1.phi(), top2.phi()));
    out_mLB->push_back((solver_->l1()+solver_->j1()).mass());
    out_mLB->push_back((solver_->l2()+solver_->j2()).mass());
    if ( jets.size() >= 4 )
    {
      int nUsedJet = 0;
      LV addJet2P4;
      for ( auto jet : jets )
      {
        if ( jet == selectedJet1 or jet == selectedJet2 ) continue;
        addJet2P4 += jet->p4();
        if ( ++nUsedJet > 2 ) break;
      }
      out_mAddJJ->push_back(addJet2P4.mass());
    }
  } while (false);

  event.put(cands);

  event.put(out_aux, "aux");
  event.put(out_mLL, "mLL");
  event.put(out_mLB, "mLB");
  event.put(out_mAddJJ, "mAddJJ");
  event.put(out_dphi, "dphi");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLLKinSolutionProducer);

