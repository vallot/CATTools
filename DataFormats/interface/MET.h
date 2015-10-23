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

    float unclusteredEnPx(int dir) const { if (dir > 0) return unclusteredEnUp_px_; else return unclusteredEnDown_px_;}
    float unclusteredEnPy(int dir) const { if (dir > 0) return unclusteredEnUp_py_; else return unclusteredEnDown_py_;}
    float unclusteredEnSumEt(int dir) const { if (dir > 0) return unclusteredEnUp_sumet_; else return unclusteredEnDown_sumet_;}
    float unclusteredEnPt(int dir) const { return hypotf(unclusteredEnPx(dir), unclusteredEnPy(dir));}
    float unclusteredEnPhi(int dir) const { return std::atan2(unclusteredEnPy(dir), unclusteredEnPx(dir));}

    void setUnclusteredEnUp (float x, float y, float s) {unclusteredEnUp_px_=x; unclusteredEnUp_py_=y; unclusteredEnUp_sumet_=s;}
    void setUnclusteredEnDown (float x, float y, float s) {unclusteredEnDown_px_=x; unclusteredEnDown_py_=y; unclusteredEnDown_sumet_=s;}
    
  private:
    float sumEt_;
    float unclusteredEnUp_px_, unclusteredEnUp_py_, unclusteredEnUp_sumet_;
    float unclusteredEnDown_px_, unclusteredEnDown_py_, unclusteredEnDown_sumet_;
    
  };
}

#endif
