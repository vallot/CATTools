#ifndef CATTools_Photon_H
#define CATTools_Photon_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

// Define typedefs for convenience
namespace cat {
  class Photon;
  typedef std::vector<Photon>              PhotonCollection;
  typedef edm::Ref<PhotonCollection>       PhotonRef;
  typedef edm::RefVector<PhotonCollection> PhotonRefVector;
}

namespace cat {

  class Photon : public Particle{
  public:
    Photon();
    Photon(const reco::LeafCandidate & aPhoton); 
    virtual ~Photon();

    double chargedHadronIso() const { return chargedHadronIso_;  }
    double puChargedHadronIso() const { return puChargedHadronIso_; }
    double neutralHadronIso() const { return neutralHadronIso_; }
    double photonIso() const { return photonIso_; }

  private:

    double chargedHadronIso_;
    double puChargedHadronIso_;
    double neutralHadronIso_;
    double photonIso_;

  };
}

#endif
