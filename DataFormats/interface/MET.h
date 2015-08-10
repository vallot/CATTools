#ifndef CATTools_MET_H
#define CATTools_MET_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"

// Define typedefs for convenience
namespace cat {
  class MET;
  typedef std::vector<MET>              METCollection;
  typedef edm::Ref<METCollection>       METRef;
  typedef edm::RefVector<METCollection> METRefVector;
}

namespace cat {

  class MET : public reco::LeafCandidate {
  public:
    typedef math::XYZTLorentzVector LorentzVector;
    
    MET() {};
    MET(float px, float py, float sumEt) { set(px, py, sumEt); };
    virtual ~MET() {};

    void set(float px, float py, float sumEt) {
      setP4(LorentzVector(px, py, 0, std::hypot(px, py)));
      sumEt_ = sumEt;
    }

    float sumEt() const { return sumEt_; }
    TLorentzVector tlv() const {return TLorentzVector(px(), py(), 0.0, pt());}

  private:
    float sumEt_;

  };
}

#endif
