#ifndef CATTools_MET_H
#define CATTools_MET_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// Define typedefs for convenience
namespace cat {
  class MET;
  typedef std::vector<MET>              METCollection;
  typedef edm::Ref<METCollection>       METRef;
  typedef edm::RefVector<METCollection> METRefVector;
}

namespace cat {

  class MET : public Particle{
  public:
    MET();
    MET(const reco::LeafCandidate & aMET); 
    virtual ~MET();

    float sumEt() const { return sumet_; }
    void setSumEt(float f){ sumet_ = f; }
    
  private:
    float sumet_;

  };
}

#endif
