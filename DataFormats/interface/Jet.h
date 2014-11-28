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
    float secvtxMass() const { return secvtxMass_ ; }
    float Lxy() const { return Lxy_ ; }
    float LxyErr() const { return LxyErr_; }

    /// \return the matched MC parton flavour (from the shower, used e.g. for b-tagging)
    int partonFlavour() const{ return partonFlavour_;}
    /// \return the pdgId of the matched MC parton from hard scattering (i.e. the closest quark or gluon of status == 3)
    int partonPdgId() const{ return partonPdgId_;}


    void setbtag_csv(float f) { csv_ = f;}
    void setSecVtxMass(float f) { secvtxMass_ = f;}
    void setLxy(float f) { Lxy_ = f;}
    void setLxyErr(float f) { LxyErr_ = f;}
    void setPartonFlavour(int i) { partonFlavour_ = i; }
    void setPartonPdgId(int i) { partonPdgId_ = i; }
    void setLooseId(bool id) { LooseId_ = id; } 
    void setPileupJetId(float f) { pileupJetId_ = f;}

    /* case CSVL: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.244; */
    /* case CSVM: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679; */
    /* case CSVT: return bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898; */

  private:

    bool LooseId_; 
    float pileupJetId_; 
    /// b tagging information
    float csv_;
    float secvtxMass_;
    float Lxy_;
    float LxyErr_;

    //parton flavour
    int partonFlavour_;
    int partonPdgId_;

  };
}

#endif
