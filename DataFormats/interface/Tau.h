#ifndef CATTools_Tau_H
#define CATTools_Tau_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

// Define typedefs for convenience
namespace cat {
  class Tau;
  typedef std::vector<Tau>              TauCollection;
  typedef edm::Ref<TauCollection>       TauRef;
  typedef edm::RefVector<TauCollection> TauRefVector;
}

namespace cat {

  class Tau : public Particle{
  public:
    Tau();
    Tau(const reco::LeafCandidate & aTau);
    virtual ~Tau();

  };
}

#endif
