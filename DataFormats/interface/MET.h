#ifndef CATTools_MET_H
#define CATTools_MET_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"

// Define typedefs for convenience
namespace cat {
  class MET;
  typedef std::vector<MET>              METCollection;
  typedef edm::Ref<METCollection>       METRef;
  typedef edm::RefVector<METCollection> METRefVector;
}

namespace cat {

  class MET : public Particle {
  public:

    MET() {};
    MET(const reco::LeafCandidate & aMet);
    MET(const reco::LeafCandidate & aMet, float sumEt);
    virtual ~MET() {};

    float sumEt() const { return sumEt_; }

  private:
    float sumEt_;

  };
}

#endif
