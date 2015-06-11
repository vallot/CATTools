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

    float ChargedMultiplicity() const { return ChargedMultiplicity_; }
    float NeutralMultiplicity() const { return NeutralMultiplicity_; }
    float ChargedHadronMultiplicity() const { return ChargedHadronMultiplicity_; }
    float NeutralHadronMultiplicity() const { return NeutralHadronMultiplicity_; }
    float MuonMultiplicity() const { return MuonMultiplicity_; }
    float ElectronMultiplicity() const { return ElectronMultiplicity_; }
    float PhotonMultiplicity() const { return PhotonMultiplicity_; }
    float HFEMMultiplicity() const { return HFEMMultiplicity_; }
    float HFHadronMultiplicity() const { return HFHadronMultiplicity_; }

    float NeutralEmEnergyFraction() const { return NeutralEmEnergyFraction_; }
    float NeutralHadronEnergyFraction() const { return NeutralHadronEnergyFraction_; }
    float ChargedEmEnergyFraction() const { return ChargedEmEnergyFraction_; }
    float ChargedHadronEnergyFraction() const { return ChargedHadronEnergyFraction_; }
    float HFEMEnergyFraction() const { return HFEMEnergyFraction_; }
    float HFHadronEnergyFraction() const { return HFHadronEnergyFraction_; }
    float ChargedMuEnergyFraction() const { return ChargedMuEnergyFraction_; }
    float MuonEnergyFraction() const { return MuonEnergyFraction_; }
    float PhotonEnergyFraction() const { return PhotonEnergyFraction_; }
    float ElectronEnergyFraction() const { return ElectronEnergyFraction_; }

    float CombinedSecondaryVertexBTag() const { return CombinedSecondaryVertexBTag_; }
    float TrackCountingHighPurBTag() const { return TrackCountingHighPurBTag_; }
    float JetProbabilityBTag() const { return JetProbabilityBTag_; }

    float L1FastJetJEC() const { return L1FastJetJEC_; }
    float L2L3ResJEC() const { return L2L3ResJEC_; }
    float L2RelJEC() const { return L2RelJEC_; }
    float L3AbsJEC() const { return L3AbsJEC_; }
//    float L5BottomJEC() const { return L5BottomJEC_; }
//    float L5ChartJEC() const { return L5CharmJEC_; }
//    float L5UDSJEC() const { return L5UDSJEC_; }
//    float L5GluonJEC() const { return L5GluonJEC_; }
//    Also will be added for L7Parton
    float EnergyRaw() const { return EnergyRaw_; }
    float PtRaw() const { return PtRaw_; }
  
 

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


    void setChargedMultiplicity(float f) { ChargedMultiplicity_ = f; }
    void setNeutralMultiplicity(float f) { NeutralMultiplicity_ = f;} 
    void setChargedHadronMultiplicity(float f) { ChargedHadronMultiplicity_ = f; }
    void setNeutralHadronMultiplicity(float f) { NeutralHadronMultiplicity_ = f; }
    void setMuonMultiplicity(float f) { MuonMultiplicity_ = f; }
    void setElectronMultiplicity(float f) { ElectronMultiplicity_ = f;}
    void setPhotonMultiplicity(float f) { PhotonMultiplicity_ = f; }
    void setHFEMMultiplicity(float f) { HFEMMultiplicity_ = f; }
    void setHFHadronMultiplicity(float f) { HFHadronMultiplicity_ = f; }

    void setNeutralEmEnergyFraction(float f) { NeutralEmEnergyFraction_ = f; }
    void setNeutralHadronEnergyFraction(float f) { NeutralHadronEnergyFraction_ = f; }
    void setChargedEmEnergyFraction(float f) { ChargedEmEnergyFraction_ = f; }
    void setChargedHadronEnergyFraction(float f) { ChargedHadronEnergyFraction_ = f; }
    void setHFEMEnergyFraction(float f) { HFEMEnergyFraction_ = f; }  
    void setHFHadronEnergyFraction(float f) { HFHadronEnergyFraction_ = f; }
    void setChargedMuEnergyFraction(float f) { ChargedMuEnergyFraction_ = f; }
    void setMuonEnergyFraction(float f) { MuonEnergyFraction_ = f; }
    void setPhotonEnergyFraction(float f) { PhotonEnergyFraction_ = f; }
    void setElectronEnergyFraction(float f) { ElectronEnergyFraction_ = f; }

    void setCombinedSecondaryVertexBTag(float f) { CombinedSecondaryVertexBTag_ = f; }
    void setTrackCountingHighPurBTag(float f) { TrackCountingHighPurBTag_ = f; }
    void setJetProbabilityBTag(float f) { JetProbabilityBTag_ = f; }

    void setL1FastJetJEC(float f) { L1FastJetJEC_ = f; }
    void setL2L3ResJEC(float f) { L2L3ResJEC_ = f; }
    void setL2RelJEC(float f) { L2RelJEC_ = f; }
    void setL3AbsJEC(float f) { L3AbsJEC_ = f; }
//    void setL5BottomJEC(float f) { L5BottomJEC_ = f; }
//    void setL5CharmJEC(float f) { L5CharmJEC_ = f; }
//    void setL5UDSJEC(float f) { L5UDSJEC_ = f; }
//    void setL5GluonJEC(float f) { L5GluonJEC_ = f; }
    void setEnergyRaw(float f) { EnergyRaw_ = f; }
    void setPtRaw(float f) { PtRaw_ = f; }



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

  private:

    float ChargedMultiplicity_;
    float NeutralMultiplicity_;
    float ChargedHadronMultiplicity_;
    float NeutralHadronMultiplicity_;
    float MuonMultiplicity_;
    float ElectronMultiplicity_;
    float PhotonMultiplicity_;
    float HFEMMultiplicity_;
    float HFHadronMultiplicity_;

    float NeutralEmEnergyFraction_;
    float NeutralHadronEnergyFraction_;
    float ChargedEmEnergyFraction_;
    float ChargedHadronEnergyFraction_;
    float HFEMEnergyFraction_;
    float HFHadronEnergyFraction_;
    float ChargedMuEnergyFraction_;
    float MuonEnergyFraction_;
    float PhotonEnergyFraction_;
    float ElectronEnergyFraction_;

    float CombinedSecondaryVertexBTag_;
    float TrackCountingHighPurBTag_;
    float JetProbabilityBTag_;

    float L1FastJetJEC_;
    float L2L3ResJEC_;
    float L2RelJEC_;
    float L3AbsJEC_;
//    float L5BottomJEC_;
//    float L5CharmJEC_;
//    float L5UDSJEC_;
//    float L5GluonJEC_;
    float EnergyRaw_;
    float PtRaw_;


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
