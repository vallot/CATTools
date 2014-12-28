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

//#include "DataFormats/JetReco/interface/GenJetCollection.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

//using namespace std;

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

  //produces<reco::CompositeRefCandidateCollection>();
}

void TTbarDileptonProducer::produce(edm::Event& event, const edm::EventSetup&) 
{
  //std::auto_ptr<reco::CompositeRefCandidateCollection> ttbarCands(new reco::CompositeRefCandidateCollection);
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
    LorentzVector nu1, nu2;
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
          nu1 = solver_->nu1();
          nu2 = solver_->nu2();
        }
      }
    }
    if ( quality <= -1e9 ) break; // failed to get solution
  } while (false);

//  ttbarCands->push_back(aGenTop);
  //event.put(ttbarCands);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTbarDileptonProducer);

