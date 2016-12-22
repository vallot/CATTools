#ifndef CATTools_Jet_H
#define CATTools_Jet_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"

#include <string>
//#include <boost/array.hpp>
#include <bitset>

// Define typedefs for convenience
namespace cat {
  class Jet;
  typedef std::vector<Jet>              JetCollection;
  typedef edm::Ref<JetCollection>       JetRef;
  typedef edm::RefVector<JetCollection> JetRefVector;
}

namespace cat {

  // pfJetProbabilityBJetTags
  const std::string BTAG_JP = "pfJetProbabilityBJetTags";
  //                         ICHEP16// 2015
  const double WP_BTAG_JPL = 0.245; // 0.275
  const double WP_BTAG_JPM = 0.515; // 0.545
  const double WP_BTAG_JPT = 0.760; // 0.790
  // pfCombinedInclusiveSecondaryVertexV2BJetTags
  const std::string BTAG_CSVv2 = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
  //                            Moriond17//ICHEP16//2015
  const double WP_BTAG_CSVv2L = 0.5426; // 0.460; // 0.605
  const double WP_BTAG_CSVv2M = 0.8484; // 0.800; // 0.89
  const double WP_BTAG_CSVv2T = 0.9535; // 0.935; // 0.97
  // pfCombinedMVAV2BJetTags
  const std::string BTAG_cMVAv2 = "pfCombinedMVAV2BJetTags";
  //                             Moriond17//ICHEP16
  const double WP_BTAG_cMVAv2L = -0.5884; // -0.715; // -
  const double WP_BTAG_cMVAv2M =  0.4432; //  0.185; // -
  const double WP_BTAG_cMVAv2T =  0.9432; //  0.875; // -
  // pfCombinedCvsLJetTags
  const std::string CTAG_CvsL = "pfCombinedCvsLJetTags";
  //                           Moriond17=ICHEP16
  const double WP_CTAG_CvsLL = -0.48;// -
  const double WP_CTAG_CvsLM = -0.1; // -
  const double WP_CTAG_CvsLT =  0.69; // -
  const std::string CTAG_CvsB = "pfCombinedCvsBJetTags";
  //                           Moriond17=ICHEP16
  const double WP_CTAG_CvsBL = -0.17; // -
  const double WP_CTAG_CvsBM = 0.08; // -
  const double WP_CTAG_CvsBT = -0.45; // -
  // DeepCSV
  const std::string BTAG_DeepCSV = "DeepCSV";
  const double WP_BTAG_DeepCSVL = 0.2219;
  const double WP_BTAG_DeepCSVM = 0.6324;
  const double WP_BTAG_DeepCSVT = 0.8958;

  class Jet : public Particle{
  public:
    Jet();
    Jet(const reco::LeafCandidate & aJet);
    virtual ~Jet();

    inline bool LooseId() const { return looseJetID(); }// temp for backward comp
    inline bool TightId() const { return tightJetID(); }// temp for backward comp
    bool looseJetID() const { return jetIdBits_[0]; }
    bool tightJetID() const { return jetIdBits_[1]; }
    bool tightLepVetoJetID() const { return jetIdBits_[2]; }

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

    bool CSVv2L(){ return (bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2L);}
    bool CSVv2M(){ return (bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2M);}
    bool CSVv2T(){ return (bDiscriminator(BTAG_CSVv2) > WP_BTAG_CSVv2T);}

    enum BTAGCSV_CUT { BTAGCSV_LOOSE=0, BTAGCSV_MEDIUM, BTAGCSV_TIGHT };
    float scaleFactorCSVv2(BTAGCSV_CUT op, int systDir) const;
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
    float smearedRes(int direction=0, int era = 0) const; // 0, +1, -1 for smeared, smearedUp, smearedDown
    float smearedResUp(int era = 0) const { return smearedRes(+1, era); };
    float smearedResDown(int era = 0) const { return smearedRes(-1, era); };

    float qgLikelihood() const { return qgLikelihood_; }
    
    void setLooseJetID(bool id) { jetIdBits_[0] = id; }
    void setTightJetID(bool id) { jetIdBits_[1] = id; }
    void setTightLepVetoJetID(bool id) { jetIdBits_[2] = id; }

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

    std::bitset<3> jetIdBits_; // loose,tight,tightLepVeto

    float pileupJetId_;
    float chargedEmEnergyFraction_;

    /// b tagging discriminators
    std::vector<std::pair<std::string, float> >  pairDiscriVector_;
    //float pairDiscriVector_;
     /// b tagging information
    float vtxMass_;
    unsigned short vtxNtracks_;
    float vtx3DVal_;
    float vtx3DSig_;

    //parton flavour
    short partonFlavour_;
    short hadronFlavour_;
    short partonPdgId_;

    float shiftedEnDown_;
    float shiftedEnUp_;

    float fJER_, fJERUp_, fJERDown_;

    float qgLikelihood_;

  };
}

#endif
