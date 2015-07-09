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

  typedef reco::PFJet::Specific PFSpecific;

  class Jet : public Particle{
  public:
    Jet();
    Jet(const reco::LeafCandidate & aJet); 
    virtual ~Jet();

    bool LooseId() const { return LooseId_; }
    bool TightId() const { return TightId_; }
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

    // ---- PF Jet specific information ----
    float photonEnergy () const {return pfSpecific().mPhotonEnergy;}
    float photonEnergyFraction () const {return photonEnergy()/correctedJet("Uncorrected").energy();}
    float electronEnergy () const {return pfSpecific().mElectronEnergy;}
    float electronEnergyFraction () const {return electronEnergy()/correctedJet("Uncorrected").energy();}
    float muonEnergy () const {return pfSpecific().mMuonEnergy;}
    float muonEnergyFraction () const {return muonEnergy()/correctedJet("Uncorrected").energy();}
    float HFHadronEnergy () const {return pfSpecific().mHFHadronEnergy;}
    float HFHadronEnergyFraction () const {return HFHadronEnergy()/correctedJet("Uncorrected").energy();}
    float HFEMEnergy () const {return pfSpecific().mHFEMEnergy;}
    float HFEMEnergyFraction () const {return HFEMEnergy()/correctedJet("Uncorrected").energy();}

    int chargedHadronMultiplicity () const {return pfSpecific().mChargedHadronMultiplicity;}
    int neutralHadronMultiplicity () const {return pfSpecific().mNeutralHadronMultiplicity;}
    int photonMultiplicity () const {return pfSpecific().mPhotonMultiplicity;}
    int electronMultiplicity () const {return pfSpecific().mElectronMultiplicity;}
    int HFHadronMultiplicity () const {return pfSpecific().mHFHadronMultiplicity;}
    int HFEMMultiplicity () const {return pfSpecific().mHFEMMultiplicity;}
    float chargedMuEnergy () const {return pfSpecific().mChargedMuEnergy;}
    float chargedMuEnergyFraction () const {return chargedMuEnergy()/correctedJet("Uncorrected").energy();}
    int neutralMultiplicity () const {return pfSpecific().mNeutralMultiplicity;}
    float hoEnergy () const {return pfSpecific().mHOEnergy;}
    float hoEnergyFraction () const {return hoEnergy()/correctedJet("Uncorrected").energy();}

    // ---- JPT or PF Jet specific information ----
    float chargedHadronEnergy() const {return pfSpecific().mChargedHadronEnergy; }
    float neutralHadronEnergy() const {return pfSpecific().mNeutralHadronEnergy; }
    float chargedEmEnergy() const {return pfSpecific().mChargedEmEnergy; }
    float neutralEmEnergy() const {return pfSpecific().mNeutralEmEnergy; }    
    int muonMultiplicity() const {return pfSpecific().mMuonMultiplicity; }
    int chargedMultiplicity() const {return pfSpecific().mChargedMultiplicity; }

    /// chargedHadronEnergyFraction (relative to uncorrected jet energy)
    float chargedHadronEnergyFraction() const {return chargedHadronEnergy()/correctedJet("Uncorrected").energy();}
    /// neutralHadronEnergyFraction (relative to uncorrected jet energy)
    float neutralHadronEnergyFraction() const {return neutralHadronEnergy()/correctedJet("Uncorrected").energy();}
    /// chargedEmEnergyFraction (relative to uncorrected jet energy)
    float chargedEmEnergyFraction()     const {return chargedEmEnergy()/correctedJet("Uncorrected").energy();}
    /// neutralEmEnergyFraction (relative to uncorrected jet energy)
    float neutralEmEnergyFraction()     const {return neutralEmEnergy()/correctedJet("Uncorrected").energy();}

    /// retrieve the pf specific part of the jet
    const PFSpecific& pfSpecific() const {
      if (specificPF_.empty()) throw cms::Exception("Type Mismatch") << "This PAT jet was not made from a PFJet.\n";
      return specificPF_[0];
    }
    void setPFSpecific(const PFSpecific& newPFSpecific) {
      specificPF_.push_back(newPFSpecific);
    }

    void setLooseId(bool id) { LooseId_ = id; }
    void setTightId(bool id) { TightId_ = id; }
    void setPileupJetId(float f) { pileupJetId_ = f;}

    void setVtxMass(float f) { vtxMass_ = f;}
    void setVtxNtracks(int f) { vtxNtracks_ = f;}
    void setVtx3DVal(float f) { vtx3DVal_ = f;}
    void setVtx3DSig(float f) { vtx3DSig_ = f;}

    void setPartonFlavour(int i) { partonFlavour_ = i; }
    void setHadronFlavour(int i) { hadronFlavour_ = i; }
    void setPartonPdgId(int i) { partonPdgId_ = i; }

    float bDiscriminator(const std::string &theLabel) const;
    /// get vector of paire labelname-disciValue
    const std::vector<std::pair<std::string, float> > & getPairDiscri() const {return pairDiscriVector_; }
    void bDiscriminatorPrint() const;
    
    void setBDiscriminators(const std::vector<std::pair<std::string, float> > & ids) { pairDiscriVector_ = ids; }
    void addBDiscriminatorPair(const std::pair<std::string, float> & thePair) {
      pairDiscriVector_.push_back(thePair);
    }

    void setShiftedEnDown(float f) { shiftedEnDown_ = f;}
    void setShiftedEnUp(float f) { shiftedEnUp_ = f;}
    void setSmearedRes(float f) { smearedRes_ = f;}
    void setSmearedResDown(float f) { smearedResDown_ = f;}
    void setSmearedResUp(float f) { smearedResUp_ = f;}

    float shiftedEnDown() {return  shiftedEnDown_;}
    float shiftedEnUp()   {return  shiftedEnUp_;}
    float smearedRes()    {return  smearedRes_;}
    float smearedResDown(){return  smearedResDown_;}
    float smearedResUp()  {return  smearedResUp_;}

    const reco::GenJet * genJet() const { return genJetFwdRef_.get();}
    void setGenJetRef(const edm::FwdRef<reco::GenJetCollection> & gj){ genJetFwdRef_ = gj;}
    //edm::FwdRef<reco::GenJetCollection> const & genJetFwdRef() const { return genJetFwdRef_; }

    void addJecFactorPair(const std::pair<std::string, float> & thePair) {
      jecFactor_.push_back(thePair);
    }
    float jecFactor(const std::string &theLabel) const;
    Jet correctedJet(const std::string &theLabel) const;

  private:

    std::vector<std::pair<std::string, float> >   jecFactor_;

    std::vector<PFSpecific>   specificPF_;

    edm::FwdRef<reco::GenJetCollection>  genJetFwdRef_;

    bool LooseId_; 
    bool TightId_;
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
