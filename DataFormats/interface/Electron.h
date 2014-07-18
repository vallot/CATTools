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

    double chargedHadronIso() const { return chargedHadronIso_;  }
    double puChargedHadronIso() const { return puChargedHadronIso_; }
    double neutralHadronIso() const { return neutralHadronIso_; }
    double photonIso() const { return photonIso_; }

    float mva() const { return mva_; }
  
  private:

    double chargedHadronIso_;
    double puChargedHadronIso_;
    double neutralHadronIso_;
    double photonIso_;

    float mva_;

  };
}

#endif
