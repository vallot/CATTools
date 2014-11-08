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

    float chargedHadronIso(float dR=0.3) const { 
      if( dR < 0.35){ return chargedHadronIso03_;}
      else { return chargedHadronIso04_;}
    }
    float puChargedHadronIso(float dR=0.3) const { 
      if( dR < 0.35) { return puChargedHadronIso03_;}
      else { return puChargedHadronIso04_;}
    }
    float neutralHadronIso(float dR=0.3) const { 
      if( dR < 0.35 ) { return neutralHadronIso03_;}
      else { return neutralHadronIso04_;}
    }
    float photonIso(float dR=0.3) const { 
      if( dR < 0.35 ) { return photonIso03_;}
      else { return photonIso04_;} 
    }

    float relIso(float dR=0.3 ) const {
      double elEffArea;
      if( dR < 0.35) elEffArea = elecEA03_; 
      else  elEffArea = elecEA04_; 
      return (chargedHadronIso(dR)+ std::max(0., neutralHadronIso(dR)+photonIso(dR)-elEffArea*rho_)) / pt_;
    }

    float mva() const { return mva_; }
    float scEta() const { return scEta_; }
    float dxy() const { return dxy_; }
    float dz() const { return dz_; }
    bool conversionVeto() const { return conversionVeto_; }
    bool chargeIDFull() const { return chargeIDFull_; }
    bool isPF() const { return isPF_; }

    void setmva(float i) { mva_ = i; }
    void setscEta(float i) { scEta_ = i; }
    void setdxy(float i) {  dxy_ = i; }
    void setdz(float i) {  dz_ = i; }
    void setconversionVeto(bool i) {  conversionVeto_ = i; }
    void setchargeIDFull(bool i) {  chargeIDFull_ = i; }
    void setisPF(bool i) {  isPF_ = i; }
    void setrho(float i) { rho_ = i; }

    void setChargedHadronIso03(float i) { chargedHadronIso03_ = i; }
    void setPUChargedHadronIso03(float i) { puChargedHadronIso03_ = i; }
    void setNeutralHadronIso03(float i) { neutralHadronIso03_ = i; }
    void setPhotonIso03(float i) { photonIso03_ = i; }
    void setEffArea03( bool runOnMC ) {
      ElectronEffectiveArea::ElectronEffectiveAreaTarget electronEATarget; 
      if ( runOnMC ) electronEATarget = ElectronEffectiveArea::kEleEAFall11MC;
      else electronEATarget = ElectronEffectiveArea::kEleEAData2012;
      elecEA03_ = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, scEta_, electronEATarget);
    }

    void setChargedHadronIso04(float i) { chargedHadronIso04_ = i; }
    void setPUChargedHadronIso04(float i) { puChargedHadronIso04_ = i; }
    void setNeutralHadronIso04(float i) { neutralHadronIso04_ = i; }
    void setPhotonIso04(float i) { photonIso04_ = i; }
    void setEffArea04( bool runOnMC ) {
      ElectronEffectiveArea::ElectronEffectiveAreaTarget electronEATarget; 
      if ( runOnMC ) electronEATarget = ElectronEffectiveArea::kEleEAFall11MC;
      else electronEATarget = ElectronEffectiveArea::kEleEAData2012;
      elecEA04_ = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, scEta_, electronEATarget);
    }

  private:

    float chargedHadronIso03_;
    float puChargedHadronIso03_;
    float neutralHadronIso03_;
    float photonIso03_;
    float elecEA03_;

    float chargedHadronIso04_;
    float puChargedHadronIso04_;
    float neutralHadronIso04_;
    float photonIso04_;
    float elecEA04_;

    float mva_;
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
