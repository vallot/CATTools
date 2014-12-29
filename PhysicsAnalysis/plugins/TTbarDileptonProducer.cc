#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "CATTools/DataFormats/interface/Particle.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
//#include "CATTools/DataFormats/interface/Photon.h"
#include "CATTools/DataFormats/interface/Jet.h"
//#include "CATTools/DataFormats/interface/Tau.h"
#include "CATTools/DataFormats/interface/MET.h"
//#include "CATTools/DataFormats/interface/GenJet.h"
//#include "CATTools/DataFormats/interface/GenTop.h"
//#include "CATTools/DataFormats/interface/MCParticle.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "CATTools/PhysicsAnalysis/interface/KinematicSolvers.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

//#include "DataFormats/JetReco/interface/GenJetCollection.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

using namespace std;

namespace cat {

class TTbarDileptonProducer : public edm::EDProducer 
{
public:
  TTbarDileptonProducer(const edm::ParameterSet& pset);
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

  solver_ = new TTDileptonSolver(); // A dummy solver
  //const auto solverName = pset.getParameter<std::string>("solver");
  //if      ( solverName == "CMSKINSOLVER" ) solver_ = DileptonSolvers::CMSKinSolver();
  //else if ( solverName == "MT2"          ) solver_ = DileptonSolvers::MT2Solver();
  //else if ( solverName == "MAOS"         ) solver_ = DileptonSolvers::MAOSSolver();

  produces<CRCandColl>();
  produces<int>("channel");
}

void TTbarDileptonProducer::produce(edm::Event& event, const edm::EventSetup&) 
{
  std::auto_ptr<CRCandColl> cands(new CRCandColl(5)); //ttbar, top1, top2, nu1, nu2
  auto candsRefProd = event.getRefBeforePut<CRCandColl>();

  std::auto_ptr<int> channel(new int);
  *channel = CH_NONE;

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
    LorentzVector lep1LVec, lep2LVec;
cout << "S1 : " << muonHandle->size() << ' ' << electronHandle->size() << endl;
    if ( muonHandle->size() + electronHandle->size() < 2 ) break; // Skip if it is not dilepton event
    // Pick leading leptons. Assume muons and electrons are sorted by pT.
    auto mu1 = muonHandle->begin(), mu2 = muonHandle->begin()+1;
    auto el1 = electronHandle->begin(), el2 = electronHandle->begin()+1;
    if ( electronHandle->empty() or mu2->pt() > mu1->pt() ) // mumu channel
    {
      *channel = CH_MUMU;
      lep1LVec = mu1->p4();
      lep2LVec = mu2->p4();
    }
    else if ( muonHandle->empty() or el2->pt() < mu1->pt() ) // elel channel
    {
      *channel = CH_ELEL;
      lep1LVec = el1->p4();
      lep2LVec = el2->p4();
    }
    else // Otherwise it is emu channel
    {
      *channel = CH_MUEL;
      lep1LVec = mu1->p4();
      lep2LVec = el1->p4();
    }

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
          selectedJet1 = jet1; selectedJet2 = jet2;
          nu1LVec = solver_->nu1();
          nu2LVec = solver_->nu2();
        }
      }
    }
cout << "S2: " << jetHandle->size() << " " << quality << endl;
    if ( quality <= -1e9 ) break; // failed to get solution

    CRCand& ttbar = cands->at(0);
    CRCand& top1 = cands->at(1);
    CRCand& top2 = cands->at(2);
    CRCand& nu1 = cands->at(3);
    CRCand& nu2 = cands->at(4);

    // Do the basic mother-daughter associations
    const auto& prodId = candsRefProd.id();
    const auto getter = candsRefProd.productGetter();
    ttbar.addDaughter(reco::CandidatePtr(prodId, 1, getter));
    ttbar.addDaughter(reco::CandidatePtr(prodId, 2, getter));

    if ( *channel == CH_MUMU )
    {
      top1.addDaughter(reco::CandidatePtr(muonHandle, mu1-muonHandle->begin()));
      top2.addDaughter(reco::CandidatePtr(muonHandle, mu2-muonHandle->begin()));
      nu1.setPdgId(14*mu1->charge());
      nu2.setPdgId(14*mu2->charge());
    }
    else if ( *channel == CH_ELEL )
    {
      top1.addDaughter(reco::CandidatePtr(electronHandle, el1-electronHandle->begin()));
      top2.addDaughter(reco::CandidatePtr(electronHandle, el2-electronHandle->begin()));
      nu1.setPdgId(14*el1->charge());
      nu2.setPdgId(14*el2->charge());
    }
    else if ( *channel == CH_MUEL )
    {
      top1.addDaughter(reco::CandidatePtr(muonHandle, mu1-muonHandle->begin()));
      top2.addDaughter(reco::CandidatePtr(electronHandle, el1-electronHandle->begin()));
      nu1.setPdgId(14*mu1->charge());
      nu2.setPdgId(14*el1->charge());
    }
cout << "S3 " << *channel << endl;

    top1.addDaughter(reco::CandidatePtr(prodId, 3, getter));
    top2.addDaughter(reco::CandidatePtr(prodId, 4, getter));
    top1.addDaughter(reco::CandidatePtr(jetHandle, selectedJet1-jetHandle->begin()));
    top2.addDaughter(reco::CandidatePtr(jetHandle, selectedJet2-jetHandle->begin()));

  } while (false);

  event.put(cands);
  event.put(channel, "channel");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTbarDileptonProducer);

