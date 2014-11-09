#ifndef CATTools_Electron_H
#define CATTools_Electron_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

// Define typedefs for convenience
namespace cat {
  class Electron;
  typedef std::vector<Electron>              ElectronCollection;
  typedef edm::Ref<ElectronCollection>       ElectronRef;
  typedef edm::RefVector<ElectronCollection> ElectronRefVector;
}

namespace cat {

  class Electron : public Particle{
  public:
    Electron();
    Electron(const reco::LeafCandidate & aElectron); 
    virtual ~Electron();

    float relIso(float dR=0.3 ) const {
      if( dR < 0.35) return relIso03_;
      else return relIso04_;
    }
    float mvaTrigV0() const { return mvaTrigV0_; }
    float scEta() const { return scEta_; }
    float dxy() const { return dxy_; }
    float dz() const { return dz_; }
    bool conversionVeto() const { return conversionVeto_; }
    bool chargeIDFull() const { return chargeIDFull_; }
    bool isPF() const { return isPF_; }

    void setrelIso(float dR, double chIso, double nhIso, double phIso, double AEff, double rhoIso, double ecalpt) {
      float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;
      if( dR < 0.35) relIso03_ = relIso; 
      else  relIso04_ = relIso; 
    }

    void setmvaTrigV0(float i) { mvaTrigV0_ = i; }
    void setscEta(float i) { scEta_ = i; }
    void setdxy(float i) {  dxy_ = i; }
    void setdz(float i) {  dz_ = i; }
    void setconversionVeto(bool i) {  conversionVeto_ = i; }
    void setchargeIDFull(bool i) {  chargeIDFull_ = i; }
    void setisPF(bool i) {  isPF_ = i; }
    void setrho(float i) { rho_ = i; }

  private:

    float relIso03_;
    float relIso04_;

    float mvaTrigV0_;
    float scEta_;
    float dxy_;
    float dz_;
    float rho_;

    bool conversionVeto_;
    bool chargeIDFull_;
    bool isPF_;
    
  };
}

#endif
