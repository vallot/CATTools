#ifndef CATTools_VertexInfo_H
#define CATTools_VertexInfo_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/Point3D.h"

namespace cat {

  class VertexInfo {
  public:
    typedef math::XYZPoint Point;
    VertexInfo();
    virtual ~VertexInfo();
    
    int numberOfVertex() const { return numberOfVertex_;}
    bool isValid() const { return validity_;}
    bool isFake() const { return isFake_;}
    Point position() const { return position_;}
    float chi2() const { return chi2_; }
    float ndof() const { return ndof_; }
    float normalizedChi2() const { return ndof_ != 0 ? chi2_ / ndof_ : chi2_ * 1e6; }
    float x() const { return position_.X(); }
    float y() const { return position_.Y(); }
    float z() const { return position_.Z(); }

    void setnumberOfVertex(int i) { numberOfVertex_ = i; }
    void setvalidity(bool i) { validity_ = i; }
    void setisFake(bool i) { isFake_ = i; }
    void setposition(const Point & i) { position_ = i; }
    void setchi2(float i) { chi2_ = i; }
    void setndof(float i) { ndof_ = i; }

  private:
    int numberOfVertex_;
    bool validity_,isFake_;
    Point position_;
    float chi2_, ndof_;
  };
}

#endif
