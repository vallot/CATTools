#include "FWCore/Framework/interface/Frameworkfwd.h"
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
#include "CommonTools/Utils/interface/PtComparator.h"

using namespace std;

namespace cat {

class TTbarDileptonProducer : public edm::stream::EDProducer<>
{
public:
  TTbarDileptonProducer(const edm::ParameterSet& pset);
  virtual ~TTbarDileptonProducer() {};
  void produce(edm::Event & event, const edm::EventSetup&) override;

private:
  typedef cat::Muon TMuon;
  typedef cat::Electron TElectron;;
  typedef cat::Jet TJet;
  typedef cat::MET TMET;
  edm::EDGetTokenT<edm::View<TMuon> > muonToken_;
  edm::EDGetTokenT<edm::View<TElectron> > electronToken_;
  edm::EDGetTokenT<edm::View<TJet> > jetToken_;
  edm::EDGetTokenT<edm::View<TMET> > metToken_;

private:
  bool isGoodMuon(const TMuon& mu)
  {
    if ( mu.pt() <= 20 or std::abs(mu.eta()) >= 2.4 ) return false;
    if ( mu.relIso(0.4) >= 0.12 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isVetoMuon(const TMuon& mu)
  {
    if ( mu.pt() <= 20 or std::abs(mu.eta()) >= 2.4 ) return false;
    if ( mu.relIso(0.4) >= 0.2 ) return false;
    if ( !mu.isLooseMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const TElectron& el)
  {
    if ( el.pt() <= 20 or std::abs(el.eta()) >= 2.4 ) return false;
    //if ( el.relIso(0.3) >= 0.11 ) return false;
    //if ( !el.electronID("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-medium") ) return false;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") ) return false;
    if ( !el.isPF() or !el.passConversionVeto() ) return false;
    const double scEta = std::abs(el.scEta());
    if ( scEta >= 1.4442 and scEta <= 1.566 ) return false;
    return true;
  }
  bool isVetoElectron(const TElectron& el)
  {
    if ( el.pt()  <= 20 or std::abs(el.eta()) >= 2.4 ) return false;
    //if ( !el.electronID("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-veto") ) return false;
    if ( !el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto") ) return false;
    const double scEta = std::abs(el.scEta());
    if ( scEta >= 1.4442 and scEta <= 1.566 ) return false;
    const double relIso03 = el.relIso(0.3);
    if ( scEta <= 1.4442 and relIso03 >= 0.1649 ) return false;
    if ( scEta >= 1.566  and relIso03 >= 0.2075 ) return false;
    return true;
  }

private:
  typedef reco::Candidate::LorentzVector LorentzVector;
  typedef reco::CompositePtrCandidate CRCand;
  typedef std::vector<CRCand> CRCandColl;
  typedef std::vector<double> doubles;
  enum CHANNEL {
    CH_NONE=0, CH_MUMU, CH_ELEL, CH_MUEL
  };
  std::unique_ptr<KinematicSolver> solver_;

};

}

using namespace cat;

TTbarDileptonProducer::TTbarDileptonProducer(const edm::ParameterSet& pset)
{
  muonToken_ = consumes<edm::View<TMuon> >(pset.getParameter<edm::InputTag>("muons"));
  electronToken_ = consumes<edm::View<TElectron> >(pset.getParameter<edm::InputTag>("electrons"));
  jetToken_ = consumes<edm::View<TJet> >(pset.getParameter<edm::InputTag>("jets"));
  metToken_ = consumes<edm::View<TMET> >(pset.getParameter<edm::InputTag>("mets"));

  auto solverName = pset.getParameter<std::string>("solver");
  std::transform(solverName.begin(), solverName.end(), solverName.begin(), ::toupper);
  if      ( solverName == "CMSKIN" ) solver_.reset(new CMSKinSolver());
  else if ( solverName == "DESYMassLoop" ) solver_.reset(new DESYMassLoopSolver());
  else if ( solverName == "DESYSmeared" ) solver_.reset(new DESYSmearedSolver());
  else if ( solverName == "MT2"    ) solver_.reset(new MT2Solver());
  else if ( solverName == "MAOS"   ) solver_.reset(new MAOSSolver());
  else if ( solverName == "NUWGT"  ) solver_.reset(new NuWeightSolver());
  else solver_.reset(new TTDileptonSolver()); // A dummy solver

  produces<CRCandColl>();
  produces<int>("channel");

  produces<doubles>("aux");
  produces<doubles>("mLL");
  produces<doubles>("mLB");
  produces<doubles>("mAddJJ");
  produces<doubles>("dphi");
  //produces<edm::RefVector<cat::Jet> >("topJets");
  //produces<edm::RefVector<cat::Jet> >("addJets");
}

void TTbarDileptonProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  std::auto_ptr<CRCandColl> cands(new CRCandColl);
  auto candsRefProd = event.getRefBeforePut<CRCandColl>();

  std::auto_ptr<int> channel(new int);
  *channel = CH_NONE;

  std::auto_ptr<doubles> out_aux(new doubles);
  std::auto_ptr<doubles> out_mLL(new doubles);
  std::auto_ptr<doubles> out_mLB(new doubles(2));
  std::auto_ptr<doubles> out_mAddJJ(new doubles);
  std::auto_ptr<doubles> out_dphi(new doubles);

  edm::Handle<edm::View<TMuon> > muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<edm::View<TElectron> > electronHandle;
  event.getByToken(electronToken_, electronHandle);

  std::vector<reco::CandidatePtr> leptons, vetoLeptons;
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    auto& p = muonHandle->at(i);
    if ( p.pt() < 20 or abs(p.eta()) >= 2.4 ) continue;

    reco::CandidatePtr muonPtr = reco::CandidatePtr(muonHandle, i);

    if ( isGoodMuon(p) ) leptons.push_back(muonPtr);
    if ( isVetoMuon(p) ) vetoLeptons.push_back(muonPtr);
  }
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    auto& p = electronHandle->at(i);
    if ( p.pt() < 20 or std::abs(p.eta()) > 2.4 ) continue;

    reco::CandidatePtr electronPtr = reco::CandidatePtr(electronHandle, i);

    if ( isGoodElectron(p) ) leptons.push_back(electronPtr);
    if ( isVetoElectron(p) ) vetoLeptons.push_back(electronPtr);
  }
  std::sort(leptons.begin(), leptons.end(), [](reco::CandidatePtr a, reco::CandidatePtr b){return a->pt() > b->pt();});
  std::sort(vetoLeptons.begin(), vetoLeptons.end(), [](reco::CandidatePtr a, reco::CandidatePtr b){return a->pt() > b->pt();});

  edm::Handle<edm::View<TJet> > jetHandle;
  event.getByToken(jetToken_, jetHandle);

  edm::Handle<edm::View<TMET> > metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& metP4 = metHandle->at(0).p4();

  do {
    // Build dileptons. We pick two highest pT muons and electrons.
    if ( leptons.size() < 2 ) break; // Skip if it is not dilepton event
    // Pick leading leptons.
    const reco::Candidate* lep1 = &*leptons.at(0);
    const reco::Candidate* lep2 = &*leptons.at(1);
    const int pdgId1 = std::abs(lep1->pdgId());
    const int pdgId2 = std::abs(lep2->pdgId());
    if ( pdgId1 == 13 and pdgId2 == 13 ) *channel = CH_MUMU;
    else if ( pdgId1 == 11 and pdgId2 == 11 ) *channel = CH_ELEL;
    else *channel = CH_MUEL;
    const LorentzVector& lep1LVec = lep1->p4();
    const LorentzVector& lep2LVec = lep2->p4();

    std::vector<reco::CandidatePtr> jets;
    for ( int i=0, n=jetHandle->size(); i<n; ++i )
    {
      const auto& jet = jetHandle->at(i);
      if ( jet.pt() < 30 or abs(jet.eta()) > 2.5 ) continue;
      if ( deltaR(jet.p4(), lep1->p4()) < 0.5 ) continue;
      if ( deltaR(jet.p4(), lep2->p4()) < 0.5 ) continue;

      jets.push_back(reco::CandidatePtr(jetHandle, i));
      //if ( jets.size() >= 2 ) break; // Just for testing.
    }

    // Run the solver with all jet combinations
    LorentzVector nu1LVec, nu2LVec;
    LorentzVector inputLVecs[5] = { metP4, lep1LVec, lep2LVec, };
    double quality = -1e9; // Default quality value
    reco::CandidatePtr selectedJet1, selectedJet2;
    for ( auto jet1 : jets )
    {
      inputLVecs[3] = jet1->p4();
      for ( auto jet2 : jets )
      {
        if ( jet1 == jet2 ) continue;
        inputLVecs[4] = jet2->p4();
        solver_->solve(inputLVecs);
        if ( solver_->quality() > quality )
        {
          quality = solver_->quality();
          selectedJet1 = jet1;
          selectedJet2 = jet2;
          nu1LVec = solver_->nu1();
          nu2LVec = solver_->nu2();
        }
      }
    }
    if ( quality <= -1e9 ) break; // failed to get solution
    // Redo the calculation with the selected ones
    inputLVecs[3] = selectedJet1->p4();
    inputLVecs[4] = selectedJet2->p4();
    solver_->solve(inputLVecs);
    std::copy(solver_->aux().begin(), solver_->aux().end(), std::back_inserter(*out_aux));

    cands->resize(7);

    CRCand& ttbar = cands->at(0);
    CRCand& top1 = cands->at(1);
    CRCand& top2 = cands->at(2);
    CRCand& w1 = cands->at(3);
    CRCand& w2 = cands->at(4);
    CRCand& nu1 = cands->at(5);
    CRCand& nu2 = cands->at(6);

    // Set four momentum
    nu1.setP4(nu1LVec);
    nu2.setP4(nu2LVec);
    w1.setP4(lep1LVec+nu1LVec);
    w2.setP4(lep2LVec+nu2LVec);
    top1.setP4(w1.p4()+selectedJet1->p4());
    top2.setP4(w2.p4()+selectedJet2->p4());
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
    auto prodId = candsRefProd.id();
    auto getter = candsRefProd.productGetter();
    ttbar.addDaughter(reco::CandidatePtr(prodId, 1, getter));
    ttbar.addDaughter(reco::CandidatePtr(prodId, 2, getter));
    top1.addDaughter(reco::CandidatePtr(prodId, 3, getter));
    top2.addDaughter(reco::CandidatePtr(prodId, 4, getter));

    w1.addDaughter(leptons.at(0));
    w2.addDaughter(leptons.at(1));
    nu1.setPdgId((pdgId1+1)*lep1->charge());
    nu2.setPdgId((pdgId2+1)*lep2->charge());

    top1.addDaughter(selectedJet1);
    top2.addDaughter(selectedJet2);
    w1.addDaughter(reco::CandidatePtr(prodId, 5, getter));
    w2.addDaughter(reco::CandidatePtr(prodId, 6, getter));

    out_mLL->push_back((lep1->p4()+lep2->p4()).mass());
    out_dphi->push_back(deltaPhi(top1.phi(), top2.phi()));
    out_mLB->push_back((lep1->p4()+selectedJet1->p4()).mass());
    out_mLB->push_back((lep2->p4()+selectedJet2->p4()).mass());
    if ( jets.size() >= 4 )
    {
      int nUsedJet = 0;
      LorentzVector addJet2P4;
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
  event.put(channel, "channel");

  event.put(out_aux, "aux");
  event.put(out_mLL, "mLL");
  event.put(out_mLB, "mLB");
  event.put(out_mAddJJ, "mAddJJ");
  event.put(out_dphi, "dphi");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTbarDileptonProducer);

