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
  // 2018 WP only
  // DeepCSV
  const std::string BTAG_DeepCSVb = "pfDeepCSVJetTags:probb";
  const std::string BTAG_DeepCSVbb = "pfDeepCSVJetTags:probbb";
  const std::string BTAG_DeepCSVc = "pfDeepCSVJetTags:probc";
  const std::string BTAG_DeepCSVcc = "pfDeepCSVJetTags:probcc";
  const std::string BTAG_DeepCSVudsg = "pfDeepCSVJetTags:probudsg";
  // DeepCSV BvsAll: pfDeepCSVJetTags:probb + pfDeepCSVJetTags:probbb
  const double WP_DeepCSVL = 0.1241;
  const double WP_DeepCSVM = 0.4184;
  const double WP_DeepCSVT = 0.7527;
  // DeepCSV CvsL: pfDeepCSVJetTags:probc/(pfDeepCSVJetTags:probc + pfDeepCSVJetTags:probudsg)
  const double WP_DeepCSV_CvsLL = 0.04;
  const double WP_DeepCSV_CvsLM = 0.137;
  const double WP_DeepCSV_CvsLT = 0.66;
  // DeepCSV CvsB: pfDeepCSVJetTags:probc/(pfDeepCSVJetTags:probc + pfDeepCSVJetTags:probb + pfDeepCSVJetTags:probbb)
  const double WP_DeepCSV_CvsBL = 0.35;
  const double WP_DeepCSV_CvsBM = 0.29;
  const double WP_DeepCSV_CvsBT = 0.10;

  // DeepJet
  const std::string BTAG_DeepJetb = "pfDeepFlavourJetTags:probb";
  const std::string BTAG_DeepJetbb = "pfDeepFlavourJetTags:probbb";
  const std::string BTAG_DeepJetlepb = "pfDeepFlavourJetTags:problepb";
  const std::string BTAG_DeepJetc = "pfDeepFlavourJetTags:probc";
  const std::string BTAG_DeepJetuds = "pfDeepFlavourJetTags:probuds";
  const std::string BTAG_DeepJetg = "pfDeepFlavourJetTags:probg";
  // DeepJet B: pfDeepFlavourJetTags:probb + pfDeepFlavourJetTags:probbb + pfDeepFlavourJetTags:problepb
  const double WP_DeepJetL = 0.0494;
  const double WP_DeepJetM = 0.2770;
  const double WP_DeepJetT = 0.7264;
  // DeepJet CvsL: pfDeepFlavourJetTags:probc/(pfDeepFlavourJetTags:probc+ pfDeepFlavourJetTags:probuds + pfDeepFlavourJetTags:probg)
  const double WP_DeepJet_CvsLL = 0.03;
  const double WP_DeepJet_CvsLM = 0.085;
  const double WP_DeepJet_CvsLT = 0.48;
  // DeepJet CvsB:  pfDeepFlavourJetTags:probc/(pfDeepFlavourJetTags:probc+ pfDeepFlavourJetTags:probb + pfDeepFlavourJetTags:probbb + pfDeepFlavourJetTags:problepb)
  const double WP_DeepJet_CvsBL = 0.4;
  const double WP_DeepJet_CvsBM = 0.29;
  const double WP_DeepJet_CvsBT = 0.05;

  class Jet : public Particle{
  public:
    Jet();
    Jet(const reco::LeafCandidate & aJet);
    virtual ~Jet();

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
    
    int printBDiscriminator() const;
    float vtxMass() const { return vtxMass_ ; }
    int vtxNtracks() const { return vtxNtracks_ ; }
    float vtx3DVal() const { return vtx3DVal_; }
    float vtx3DSig() const { return vtx3DSig_; }

    float shiftedEnDown() const {return  shiftedEnDown_;}
    float shiftedEnUp() const  {return  shiftedEnUp_;}
    float smearedRes(int direction=0, int era = 0) const; // 0, +1, -1 for smeared, smearedUp, smearedDown
    float smearedResUp(int era = 0) const { return smearedRes(+1, era); };
    float smearedResDown(int era = 0) const { return smearedRes(-1, era); };

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

    void setJER(float fJER, float fJERUp, float fJERDown) {
      fJER_ = fJER; fJERUp_ = fJERUp; fJERDown_ = fJERDown;
    }

    void setQGLikelihood(float f) { qgLikelihood_ = f; }

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

    float fJER_, fJERUp_, fJERDown_;

    float qgLikelihood_;

  };
}

#endif
