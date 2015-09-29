#ifndef CATTools_SecVertex_H
#define CATTools_SecVertex_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "CATTools/DataFormats/interface/Particle.h"
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
    int ipos() const { return ipos_;}
    int ineg() const { return ineg_;}

    void setLxy(float i) { lxy_ = i; }
    void setL3D(float i) { l3D_ = i; }
    void setVProb(float i) { vProb_ = i; }
    void setInts(int i, int j) { ipos_ = i; ineg_ = j;}

  private:
    float lxy_, l3D_, vProb_;
    int ipos_, ineg_;
  };
}

#endif
