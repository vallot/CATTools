#ifndef CATTools_Electron_H
#define CATTools_Electron_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// Define typedefs for convenience
namespace cat {
  class Electron;
  typedef std::vector<Electron>              ElectronCollection;
  typedef edm::Ref<ElectronCollection>       ElectronRef;
  typedef edm::RefVector<ElectronCollection> ElectronRefVector;
}

namespace cat {

  class Electron : public Particle{
  public:
    Electron();
    Electron(const reco::LeafCandidate & aElectron); 
    virtual ~Electron();

    float chargedHadronIso() const { return chargedHadronIso_;  }
    float puChargedHadronIso() const { return puChargedHadronIso_; }
    float neutralHadronIso() const { return neutralHadronIso_; }
    float photonIso() const { return photonIso_; }

    float mva() const { return mva_; }
  
  private:

    float chargedHadronIso_;
    float puChargedHadronIso_;
    float neutralHadronIso_;
    float photonIso_;

    float mva_;

  };
}

#endif
