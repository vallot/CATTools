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

    float chargedHadronIso() const { return chargedHadronIso_;  }
    float puChargedHadronIso() const { return puChargedHadronIso_; }
    float neutralHadronIso() const { return neutralHadronIso_; }
    float photonIso() const { return photonIso_; }

  private:

    float chargedHadronIso_;
    float puChargedHadronIso_;
    float neutralHadronIso_;
    float photonIso_;

  };
}

#endif
