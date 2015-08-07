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

  class MET {
  public:
    typedef math::XYZTLorentzVector LorentzVector;
    
    MET();
    MET(float px, float py, float sumEt);
    virtual ~MET();

    void set(float px, float py, float sumEt) { px_ = px; py_ = py; sumEt_ = sumEt;}

    float sumEt() const { return sumEt_; }
    float px() const { return px_; }
    float py() const { return py_; }
    float pt() const { return sqrt(px_*px_ + py_*py_); }
    float phi() const { return px_ == 0.0 && py_ == 0.0 ? 0.0 : TMath::ATan2(py_,px_); }
    TLorentzVector tlv() const {return TLorentzVector(px_, py_, 0.0, this->pt() );}
    LorentzVector p4()  const {return LorentzVector(px_, py_, 0.0, this->pt() ) ;}
  private:

    float sumEt_, px_, py_;

  };
}

#endif
