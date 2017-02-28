#ifndef CATTools_MET_H
#define CATTools_MET_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"

// Define typedefs for convenience
namespace cat {
  class MET;
  typedef std::vector<MET>              METCollection;
  typedef edm::Ref<METCollection>       METRef;
  typedef edm::RefVector<METCollection> METRefVector;
  
}

namespace cat {

  class MET : public Particle {
  public:

    MET() {};
    MET(const reco::LeafCandidate & aMet);
    MET(const reco::LeafCandidate & aMet, float sumEt);
    virtual ~MET() {};

    float sumEt() const { return sumEt_; }
    float rawMET() const {return rawMET_;}
    float unclusteredEnPx(int dir) const { if (dir == 1) return unclusteredEnUp_px_; else if (dir == -1) return unclusteredEnDown_px_; else return 0.;}
    float unclusteredEnPy(int dir) const { if (dir == 1) return unclusteredEnUp_py_; else if (dir == -1) return unclusteredEnDown_py_; else return 0.;}
    float unclusteredEnSumEt(int dir) const { if (dir == 1) return unclusteredEnUp_sumet_; else if (dir == -1) return unclusteredEnDown_sumet_; else return 0.;}
    float unclusteredEnPt(int dir) const { return hypotf(unclusteredEnPx(dir), unclusteredEnPy(dir));}
    float unclusteredEnPhi(int dir) const { return std::atan2(unclusteredEnPy(dir), unclusteredEnPx(dir));}

    float JetEnPx(int dir) const { if (dir == 1) return jetEnUp_px_; else if (dir == -1) return jetEnDown_px_;  else return 0.;}
    float JetEnPy(int dir) const { if (dir == 1) return jetEnUp_py_; else if (dir == -1) return jetEnDown_py_; else return 0.;}
    float JetEnSumEt(int dir) const { if (dir == 1) return jetEnUp_sumet_; else if (dir == -1) return jetEnDown_sumet_; else return 0.;}
    float JetResPx(int dir) const { if (dir == 1) return jetResUp_px_; else if (dir == -1) return jetResDown_px_; else return 0.;}
    float JetResPy(int dir) const { if (dir == 1) return jetResUp_py_; else if (dir == -1) return jetResDown_py_; else return 0.;}
    float JetResSumEt(int dir) const { if (dir == 1) return jetResUp_sumet_; else if (dir == -1) return jetResDown_sumet_; else return 0.;}
    
    float XYShiftedPx() const {return xyshift_x;}
    float XYShiftedPy() const {return xyshift_y;}
    float XYShiftedsumEt() const {return xyshift_sumet;}

    void setRawMET(float m) {rawMET_=m;}
    void setUnclusteredEnUp (float x, float y, float s) {unclusteredEnUp_px_=x; unclusteredEnUp_py_=y; unclusteredEnUp_sumet_=s;}
    void setUnclusteredEnDown (float x, float y, float s) {unclusteredEnDown_px_=x; unclusteredEnDown_py_=y; unclusteredEnDown_sumet_=s;}
    void setJetEnUp (float x, float y, float s) {jetEnUp_px_=x; jetEnUp_py_=y; jetEnUp_sumet_=s;}
    void setJetEnDown (float x, float y, float s) {jetEnDown_px_=x; jetEnDown_py_=y; jetEnDown_sumet_=s;}
    void setJetResUp (float x, float y, float s) {jetResUp_px_=x; jetResUp_py_=y; jetResUp_sumet_=s;}
    void setJetResDown (float x, float y, float s) {jetResDown_px_=x; jetResDown_py_=y; jetResDown_sumet_=s;}
    void setXYShiftedMET (float x, float y, float s){ xyshift_x=x, xyshift_y=y, xyshift_sumet=s;}
  private:
    float sumEt_;
    float rawMET_;
    float unclusteredEnUp_px_, unclusteredEnUp_py_, unclusteredEnUp_sumet_;
    float unclusteredEnDown_px_, unclusteredEnDown_py_, unclusteredEnDown_sumet_;
    float jetEnUp_px_, jetEnUp_py_, jetEnUp_sumet_;
    float jetEnDown_px_, jetEnDown_py_, jetEnDown_sumet_;
    float jetResUp_px_, jetResUp_py_, jetResUp_sumet_;
    float jetResDown_px_, jetResDown_py_, jetResDown_sumet_;
    float xyshift_x, xyshift_y, xyshift_sumet;

  };
}

#endif
