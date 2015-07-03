#ifndef CATTools_Jet_H
#define CATTools_Jet_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <string>
#include <boost/array.hpp>

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
    float vtxMass() const { return vtxMass_ ; }
    int vtxNtracks() const { return vtxNtracks_ ; }
    float vtx3DVal() const { return vtx3DVal_; }
    float vtx3DSig() const { return vtx3DSig_; }

    // return the matched MC parton flavour (from the shower, used e.g. for b-tagging)
    int partonFlavour() const{ return partonFlavour_;}
    int hadronFlavour() const{ return hadronFlavour_;}
    // pdgId of the matched MC parton from hard scattering (i.e. the closest quark or gluon of status == 3)
    int partonPdgId() const{ return partonPdgId_;}


    void setLooseId(bool id) { LooseId_ = id; } 
    void setPileupJetId(float f) { pileupJetId_ = f;}

    void setVtxMass(float f) { vtxMass_ = f;}
    void setVtxNtracks(int f) { vtxNtracks_ = f;}
    void setVtx3DVal(float f) { vtx3DVal_ = f;}
    void setVtx3DSig(float f) { vtx3DSig_ = f;}

    void setPartonFlavour(int i) { partonFlavour_ = i; }
    void setHadronFlavour(int i) { hadronFlavour_ = i; }
    void setPartonPdgId(int i) { partonPdgId_ = i; }

    float bDiscriminator(const std::string &theLabel) const;
    void setBDiscriminators(const std::vector<std::pair<std::string, float> > & ids) { pairDiscriVector_ = ids; }
    void addBDiscriminatorPair(const std::pair<std::string, float> & thePair) {
      pairDiscriVector_.push_back(thePair);
    }

    void setShiftedEnDown(float f) { shiftedEnDown_ = f;}
    void setShiftedEnUp(float f) { shiftedEnUp_ = f;}
    void setSmearedRes(float f) { smearedRes_ = f;}
    void setSmearedResDown(float f) { smearedResDown_ = f;}
    void setSmearedResUp(float f) { smearedResUp_ = f;}

    float shiftedEnDown() const {return  shiftedEnDown_;}
    float shiftedEnUp() const  {return  shiftedEnUp_;}
    float smearedRes() const {return  smearedRes_;}
    float smearedResDown() const {return  smearedResDown_;}
    float smearedResUp() const {return  smearedResUp_;}

    const reco::GenJet * genJet() const { return genJetFwdRef_.get();}
    void setGenJetRef(const edm::FwdRef<reco::GenJetCollection> & gj){ genJetFwdRef_ = gj;}
    //edm::FwdRef<reco::GenJetCollection> const & genJetFwdRef() const { return genJetFwdRef_; }

  private:
    edm::FwdRef<reco::GenJetCollection>  genJetFwdRef_;

    bool LooseId_; 
    float pileupJetId_;

    /// b tagging discriminators
    std::vector<std::pair<std::string, float> >  pairDiscriVector_;
     /// b tagging information
    float vtxMass_;
    int vtxNtracks_;
    float vtx3DVal_;
    float vtx3DSig_;

    //parton flavour
    int partonFlavour_;
    int hadronFlavour_;
    int partonPdgId_;

    float shiftedEnDown_;
    float shiftedEnUp_;
    float smearedRes_;
    float smearedResDown_;
    float smearedResUp_;

  };
}

#endif
