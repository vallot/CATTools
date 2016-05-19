#ifndef CATTools_SecVertex_H
#define CATTools_SecVertex_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

// Define typedefs for convenience
namespace cat {
  class SecVertex;
  typedef std::vector<SecVertex>              SecVertexCollection;
  typedef edm::Ref<SecVertexCollection>       SecVertexRef;
  typedef edm::RefVector<SecVertexCollection> SecVertexRefVector;
}

namespace cat {

  class SecVertex : public reco::VertexCompositeCandidate{
  public:
    SecVertex();
    SecVertex(reco::VertexCompositeCandidate & aSecVertex);
    virtual ~SecVertex();

    float lxy() const { return lxy_;}
    float l3D() const { return l3D_;}
    float vProb() const { return vProb_;}
    int leptonID1() const { return leptonID1_;}
    int trackQuality1() const { return trackQuality1_;}
    int leptonID2() const { return leptonID2_;}
    int trackQuality2() const { return trackQuality2_;}

    void setLxy(float i) { lxy_ = i; }
    void setL3D(float i) { l3D_ = i; }
    void setVProb(float i) { vProb_ = i; }
    void setLeptonID(int i, int j=-1) { leptonID1_= i; leptonID2_ = j;}
    void setTrackQuality(int i, int j=-1) { trackQuality1_= i; trackQuality2_ = j;}
    void setMCMatch(bool flag) { isMCMatch_ = flag ; }
    void setJetDR( float dR ) { jetDR_ = dR; }
    void setLegDR( float dR ) { legDR_ = dR; }

    float dca() const { return dca_;}// distance of closest approach
    float dca(int i) const { 
      if ( i==0) return dca_;
      else if ( i==1 ) return dca2_;
      else if ( i==2 ) return dca3_;
      else return -9;
    }
    void set_dca(float i) { dca_ = i; }
    void set_dca(int i, float dca) { 
      if ( i==0 ) dca_ = dca; 
      else if (i==1) dca2_= dca ; 
      else if (i==2) dca3_ = dca;
    }

    float cxPtHypot() const { return cxPtHypot_;}// crossing point hypot
    void set_cxPtHypot(float i) { cxPtHypot_ = i; }
    float cxPtAbs() const { return cxPtAbs_;}// crossing point abs
    void set_cxPtAbs(float i) { cxPtAbs_ = i; }
  private:
    float lxy_, l3D_, vProb_, dca_,dca2_,dca3_, cxPtHypot_, cxPtAbs_, jetDR_, legDR_ ;
    //int ipos_, ineg_;
    int leptonID1_, trackQuality1_;
    int leptonID2_, trackQuality2_;
    bool isMCMatch_;

  };
}

#endif
