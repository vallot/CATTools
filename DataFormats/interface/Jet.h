#ifndef CATTools_Jet_H
#define CATTools_Jet_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// Define typedefs for convenience
namespace cat {
  class Jet;
  typedef std::vector<Jet>              JetCollection;
  typedef edm::Ref<JetCollection>       JetRef;
  typedef edm::RefVector<JetCollection> JetRefVector;
}

namespace cat {

  class Jet : public Particle{
  public:
    Jet();
    Jet(const reco::LeafCandidate & aJet); 
    virtual ~Jet();

    bool LooseId() const { return LooseId_; }
    float pileupJetId() const { return pileupJetId_; }

    /// \return secondary vertex b-tagging information
    // combinedSecondaryVertexBJetTags
    float btag_csv() const{ return csv_;}
    float vtxMass() const { return vtxMass_ ; }
    int vtxNtracks() const { return vtxNtracks_ ; }
    float vtx3DVal() const { return vtx3DVal_; }
    float vtx3DSig() const { return vtx3DSig_; }

    // return the matched MC parton flavour (from the shower, used e.g. for b-tagging)
    int partonFlavour() const{ return partonFlavour_;}
    // pdgId of the matched MC parton from hard scattering (i.e. the closest quark or gluon of status == 3)
    int partonPdgId() const{ return partonPdgId_;}


    void setLooseId(bool id) { LooseId_ = id; } 
    void setPileupJetId(float f) { pileupJetId_ = f;}

    void setBtag_csv(float f) { csv_ = f;}
    void setVtxMass(float f) { vtxMass_ = f;}
    void setVtxNtracks(int f) { vtxNtracks_ = f;}
    void setVtx3DVal(float f) { vtx3DVal_ = f;}
    void setVtx3DSig(float f) { vtx3DSig_ = f;}

    void setPartonFlavour(int i) { partonFlavour_ = i; }
    void setPartonPdgId(int i) { partonPdgId_ = i; }

    /* case CSVL: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.244; */
    /* case CSVM: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679; */
    /* case CSVT: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898; */

  private:

    bool LooseId_; 
    float pileupJetId_;
     /// b tagging information
    float csv_;
    float vtxMass_;
    int vtxNtracks_;
    float vtx3DVal_;
    float vtx3DSig_;

    //parton flavour
    int partonFlavour_;
    int partonPdgId_;

  };
}

#endif
