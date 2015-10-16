#ifndef CATTools_Jet_H
#define CATTools_Jet_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

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

    bool LooseId() const { return looseJetID_; }// temp for backward comp
    bool TightId() const { return tightJetID_; }// temp for backward comp
    bool looseJetID() const { return looseJetID_; }
    bool tightJetID() const { return tightJetID_; }
    bool tightLepVetoJetID() const { return tightLepVetoJetID_; }

    float pileupJetId() const { return pileupJetId_; }
    float chargedEmEnergyFraction() const { return chargedEmEnergyFraction_; }

    // return the matched MC parton flavour (from the shower, used e.g. for b-tagging)
    int partonFlavour() const{ return partonFlavour_;}
    int hadronFlavour() const{ return hadronFlavour_;}
    // pdgId of the matched MC parton from hard scattering (i.e. the closest quark or gluon of status == 3)
    int partonPdgId() const{ return partonPdgId_;}
    const reco::GenJet * genJet() const { return genJetFwdRef_.get();}

    float bDiscriminator(const std::string &theLabel) const;
    const std::vector<std::pair<std::string, float> > & getPairDiscri() const {return pairDiscriVector_; }

    bool CSVv2L(){ return (bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.605);}
    bool CSVv2M(){ return (bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.890);}
    bool CSVv2T(){ return (bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.970);}

    enum BTAGCSV_CUT { BTAGCSV_LOOSE=0, BTAGCSV_MEDIUM, BTAGCSV_TIGHT };
    enum JETFLAV { JETFLAV_B=0, JETFLAV_C=1, JETFLAV_LIGHT=2 };
    float scaleFactorCSVv2(BTAGCSV_CUT op, int systDir, JETFLAV flav) const;
    // op 0 = loose, 1 = medium, 2 = tight
    // sys 0 = central, 1 = up, -1 = down
    // flav 0 = b, 1 = c, 2 = light
    
    void bDiscriminatorPrint() const;
    // combinedSecondaryVertexBJetTags
    float vtxMass() const { return vtxMass_ ; }
    int vtxNtracks() const { return vtxNtracks_ ; }
    float vtx3DVal() const { return vtx3DVal_; }
    float vtx3DSig() const { return vtx3DSig_; }

    float shiftedEnDown() const {return  shiftedEnDown_;}
    float shiftedEnUp() const  {return  shiftedEnUp_;}
    float smearedRes(int direction=0) const; // 0, +1, -1 for smeared, smearedUp, smearedDown
    float smearedResUp() const { return smearedRes(+1); };
    float smearedResDown() const { return smearedRes(-1); };

    
    void setLooseJetID(bool id) { looseJetID_ = id; }
    void setTightJetID(bool id) { tightJetID_ = id; }
    void setTightLepVetoJetID(bool id) { tightLepVetoJetID_ = id; }

    void setPileupJetId(float f) { pileupJetId_ = f;}
    void setChargedEmEnergyFraction(float f) { chargedEmEnergyFraction_ = f; }

    void setVtxMass(float f) { vtxMass_ = f;}
    void setVtxNtracks(int f) { vtxNtracks_ = f;}
    void setVtx3DVal(float f) { vtx3DVal_ = f;}
    void setVtx3DSig(float f) { vtx3DSig_ = f;}

    void setPartonFlavour(int i) { partonFlavour_ = i; }
    void setHadronFlavour(int i) { hadronFlavour_ = i; }
    void setPartonPdgId(int i) { partonPdgId_ = i; }
    void setGenJetRef(const edm::FwdRef<reco::GenJetCollection> & gj){ genJetFwdRef_ = gj;}

    void setBDiscriminators(const std::vector<std::pair<std::string, float> > & ids) { pairDiscriVector_ = ids; }
    void addBDiscriminatorPair(const std::pair<std::string, float> & thePair) {pairDiscriVector_.push_back(thePair);}
    /* void addBDiscriminatorPair(float f) { pairDiscriVector_ = f;} */
    /* float bDiscriminator() const {return  pairDiscriVector_;} */

    void setShiftedEnDown(float f) { shiftedEnDown_ = f;}
    void setShiftedEnUp(float f) { shiftedEnUp_ = f;}

  private:

    edm::FwdRef<reco::GenJetCollection>  genJetFwdRef_;

    bool looseJetID_;
    bool tightJetID_;
    bool tightLepVetoJetID_;

    float pileupJetId_;
    float chargedEmEnergyFraction_;

    /// b tagging discriminators
    std::vector<std::pair<std::string, float> >  pairDiscriVector_;
    //float pairDiscriVector_;
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

  };
}

#endif
