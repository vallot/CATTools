#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "CATTools/PhysicsAnalysis/interface/KinematicSolvers.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

using namespace std;

namespace cat {

class TTbarDileptonProducer : public edm::EDProducer 
{
public:
  TTbarDileptonProducer(const edm::ParameterSet& pset);
  virtual ~TTbarDileptonProducer();
  void produce(edm::Event & event, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<cat::MuonCollection> muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<cat::METCollection> metToken_;

private:
  typedef reco::Candidate::LorentzVector LorentzVector;
  typedef reco::CompositePtrCandidate CRCand;
  typedef std::vector<CRCand> CRCandColl;
  typedef std::vector<double> doubles;
  enum CHANNEL {
    CH_NONE=0, CH_MUMU, CH_ELEL, CH_MUEL
  };
  KinematicSolver* solver_;

};

}

using namespace cat;

TTbarDileptonProducer::TTbarDileptonProducer(const edm::ParameterSet& pset)
{
  muonToken_ = consumes<cat::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));
  electronToken_ = consumes<cat::ElectronCollection>(pset.getParameter<edm::InputTag>("electrons"));
  jetToken_ = consumes<cat::JetCollection>(pset.getParameter<edm::InputTag>("jets"));
  metToken_ = consumes<cat::METCollection>(pset.getParameter<edm::InputTag>("mets"));

  //pseudoTopToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InpuTag>("pseudoTop"));

  auto solverName = pset.getParameter<std::string>("solver");
  std::transform(solverName.begin(), solverName.end(), solverName.begin(), ::toupper);
  if      ( solverName == "CMSKIN" ) solver_ = new CMSKinSolver();
  else if ( solverName == "MT2"    ) solver_ = new MT2Solver();
  else if ( solverName == "MAOS"   ) solver_ = new MAOSSolver();
  else solver_ = new TTDileptonSolver(); // A dummy solver

  produces<CRCandColl>();
  produces<int>("channel");

  produces<doubles>("mLL");
  produces<doubles>("mLB");
  produces<doubles>("mAddJJ");
  produces<doubles>("dphi");
  //produces<edm::RefVector<cat::Jet> >("topJets");
  //produces<edm::RefVector<cat::Jet> >("addJets");
}

TTbarDileptonProducer::~TTbarDileptonProducer()
{
  if ( solver_ ) delete solver_;
}

void TTbarDileptonProducer::produce(edm::Event& event, const edm::EventSetup&) 
{
  std::auto_ptr<CRCandColl> cands(new CRCandColl);
  auto candsRefProd = event.getRefBeforePut<CRCandColl>();

  std::auto_ptr<int> channel(new int);
  *channel = CH_NONE;

  std::auto_ptr<doubles> out_mLL(new doubles);
  std::auto_ptr<doubles> out_mLB(new doubles(2));
  std::auto_ptr<doubles> out_mAddJJ(new doubles);
  std::auto_ptr<doubles> out_dphi(new doubles);

  edm::Handle<cat::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<cat::ElectronCollection> electronHandle;
  event.getByToken(electronToken_, electronHandle);

  edm::Handle<cat::JetCollection> jetHandle;
  event.getByToken(jetToken_, jetHandle);

  edm::Handle<cat::METCollection> metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& met = metHandle->at(0);

  do {
    // Build dileptons. We pick two highest pT muons and electrons.
    const reco::Candidate* lep1 = 0, * lep2 = 0;
    const int nMuon = muonHandle->size();
    const int nElectron = electronHandle->size();
    if ( nMuon + nElectron < 2 ) break; // Skip if it is not dilepton event
    // Pick leading leptons. Assume muons and electrons are sorted by pT.
    auto mu1 = muonHandle->begin(), mu2 = muonHandle->begin()+1;
    auto el1 = electronHandle->begin(), el2 = electronHandle->begin()+1;
    if ( electronHandle->empty() or (nMuon>=2 and mu2->pt() > el1->pt()) ) // mumu channel
    {
      *channel = CH_MUMU;
      lep1 = &*mu1; lep2 = &*mu2;
    }
    else if ( muonHandle->empty() or (nElectron>=2 and el2->pt() > mu1->pt()) ) // elel channel
    {
      *channel = CH_ELEL;
      lep1 = &*el1; lep2 = &*el2;
    }
    else // Otherwise it is emu channel
    {
      *channel = CH_MUEL;
      lep1 = &*mu1; lep2 = &*el1;
    }
    const LorentzVector& lep1LVec = lep1->p4();
    const LorentzVector& lep2LVec = lep2->p4();

    // Run the solver with all jet combinations
    LorentzVector nu1LVec, nu2LVec;
    LorentzVector inputLVecs[5] = { lep1LVec, lep2LVec, met.p4(), };
    double quality = -1e9; // Default quality value
    auto selectedJet1 = jetHandle->end(), selectedJet2 = jetHandle->end();
    for ( auto jet1 = jetHandle->begin(); jet1 != jetHandle->end(); ++jet1 )
    {
      inputLVecs[3] = jet1->p4();
      for ( auto jet2 = jetHandle->begin(); jet2 != jetHandle->end(); ++jet2 )
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

    if ( *channel == CH_MUMU )
    {
      w1.addDaughter(reco::CandidatePtr(muonHandle, mu1-muonHandle->begin()));
      w2.addDaughter(reco::CandidatePtr(muonHandle, mu2-muonHandle->begin()));
      nu1.setPdgId(14*mu1->charge());
      nu2.setPdgId(14*mu2->charge());
    }
    else if ( *channel == CH_ELEL )
    {
      w1.addDaughter(reco::CandidatePtr(electronHandle, el1-electronHandle->begin()));
      w2.addDaughter(reco::CandidatePtr(electronHandle, el2-electronHandle->begin()));
      nu1.setPdgId(12*el1->charge());
      nu2.setPdgId(12*el2->charge());
    }
    else if ( *channel == CH_MUEL )
    {
      w1.addDaughter(reco::CandidatePtr(muonHandle, mu1-muonHandle->begin()));
      w2.addDaughter(reco::CandidatePtr(electronHandle, el1-electronHandle->begin()));
      nu1.setPdgId(14*mu1->charge());
      nu2.setPdgId(12*el1->charge());
    }

    top1.addDaughter(reco::CandidatePtr(jetHandle, selectedJet1-jetHandle->begin()));
    top2.addDaughter(reco::CandidatePtr(jetHandle, selectedJet2-jetHandle->begin()));
    w1.addDaughter(reco::CandidatePtr(prodId, 5, getter));
    w2.addDaughter(reco::CandidatePtr(prodId, 6, getter));

    out_mLL->push_back((lep1->p4()+lep2->p4()).mass());
    out_dphi->push_back(deltaPhi(top1.phi(), top2.phi()));
    out_mLB->push_back((lep1->p4()+selectedJet1->p4()).mass());
    out_mLB->push_back((lep2->p4()+selectedJet2->p4()).mass());
    if ( jetHandle->size() >= 4 )
    {
      int nUsedJet = 0;
      LorentzVector addJet2P4;
      for ( auto jet = jetHandle->begin(); jet != jetHandle->end(); ++jet )
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

  event.put(out_mLL, "mLL");
  event.put(out_mLB, "mLB");
  event.put(out_mAddJJ, "mAddJJ");
  event.put(out_dphi, "dphi");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTbarDileptonProducer);

