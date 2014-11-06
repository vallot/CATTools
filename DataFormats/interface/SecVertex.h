#ifndef CATTools_SecVertex_H
#define CATTools_SecVertex_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// Define typedefs for convenience
namespace cat {
  class SecVertex;
  typedef std::vector<SecVertex>              SecVertexCollection;
  typedef edm::Ref<SecVertexCollection>       SecVertexRef;
  typedef edm::RefVector<SecVertexCollection> SecVertexRefVector;
}

namespace cat {

  class SecVertex : public Particle{
  public:
    SecVertex();
    SecVertex(const reco::LeafCandidate & aSecVertex); 
    virtual ~SecVertex();
    
    double lxy() const { return lxy_;}
    double l3D() const { return l3D_;}
    double vProb() const { return vProb_;}

    void setLxy(float i) { lxy_ = i; }
    void setL3D(float i) { l3D_ = i; }
    void setVProb(float i) { vProb_ = i; }

  private:
    float lxy_, l3D_, vProb_;
    int int_pos, int_neg;
  };
}

#endif
