#ifndef CATTools_Electron_H
#define CATTools_Electron_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

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
    float scEta() const { return scEta_; }
    float dxy() const { return dxy_; }
    float dz() const { return dz_; }
    bool passConversionVeto() const { return passConversionVeto_; }
    bool isGsfCtfScPixChargeConsistent() const { return isGsfCtfScPixChargeConsistent_; }
    bool isPF() const { return isPF_; }

    float chargedHadronIso(float dR=0.3) const {
      if( dR == 0.3){ return chargedHadronIso03_;}
      else if( dR == 0.4) { return chargedHadronIso04_;}
      else return -1.0;
    }
    float puChargedHadronIso(float dR=0.3) const {
      if( dR == 0.3) { return puChargedHadronIso03_;}
      else if( dR == 0.4) { return puChargedHadronIso04_;}
      else return -1.0;
    }
    float neutralHadronIso(float dR=0.3) const {
      if( dR == 0.3) { return neutralHadronIso03_;}
      else if( dR == 0.4) { return neutralHadronIso04_;}
      else return -1.0;
    }
    float photonIso(float dR=0.3) const {
      if( dR == 0.3) { return photonIso03_;}
      else if( dR == 0.4) { return photonIso04_;}
      else return -1.0;
    }

    float absIso(float dR=0.3, float dBetaFactor=0.5) const{

      if(dBetaFactor>0 && puChargedHadronIso(dR)<0) return -1;

      float neutralIso = neutralHadronIso(dR) + photonIso(dR);
      float corNeutralIso = neutralIso - dBetaFactor * puChargedHadronIso(dR);

      float charged = chargedHadronIso(dR);

      return charged + ( corNeutralIso>0 ? corNeutralIso : 0 ) ;
    }

    bool mcMatched() const { return mcMatched_; }

    void setrelIso(double dR, double chIso, double nhIso, double phIso, double AEff, double rhoIso, double ecalpt)
    {
      float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;
      if( dR < 0.35) relIso03_ = relIso; 
      else  relIso04_ = relIso; 
    }
    void setscEta(float i) { scEta_ = i; }
    void setdxy(float i) {  dxy_ = i; }
    void setdz(float i) {  dz_ = i; }
    void setPassConversionVeto(bool i) {  passConversionVeto_ = i; }
    void setIsGsfCtfScPixChargeConsistent(bool i) {  isGsfCtfScPixChargeConsistent_ = i; }
    void setisPF(bool i) {  isPF_ = i; }
    void setrho(float i) { rho_ = i; }

    void setChargedHadronIso03(float i) { chargedHadronIso03_ = i; }
    void setPUChargedHadronIso03(float i) { puChargedHadronIso03_ = i; }
    void setNeutralHadronIso03(float i) { neutralHadronIso03_ = i; }
    void setPhotonIso03(float i) { photonIso03_ = i; }

    void setChargedHadronIso04(float i) { chargedHadronIso04_ = i; }
    void setPUChargedHadronIso04(float i) { puChargedHadronIso04_ = i; }
    void setNeutralHadronIso04(float i) { neutralHadronIso04_ = i; }
    void setPhotonIso04(float i) { photonIso04_ = i; }

    void setMCMatched(bool m) { mcMatched_ = m; }

    float cutBasedElectronIDveto() const { return cutBasedElectronIDveto_; }
    void setcutBasedElectronIDveto(float i) { cutBasedElectronIDveto_ = i; }
    float cutBasedElectronIDloose() const { return cutBasedElectronIDloose_; }
    void setcutBasedElectronIDloose(float i) { cutBasedElectronIDloose_ = i; }
    float cutBasedElectronIDmedium() const { return cutBasedElectronIDmedium_; }
    void setcutBasedElectronIDmedium(float i) { cutBasedElectronIDmedium_ = i; }
    float cutBasedElectronIDtight() const { return cutBasedElectronIDtight_; }
    void setcutBasedElectronIDtight(float i) { cutBasedElectronIDtight_ = i; }

  private:

    float relIso03_;
    float relIso04_;

    float cutBasedElectronIDveto_;
    float cutBasedElectronIDloose_;
    float cutBasedElectronIDmedium_;
    float cutBasedElectronIDtight_;

    float scEta_;
    float dxy_;
    float dz_;
    float rho_;
    float chargedHadronIso03_;
    float puChargedHadronIso03_;
    float neutralHadronIso03_;
    float photonIso03_;

    float chargedHadronIso04_;
    float puChargedHadronIso04_;
    float neutralHadronIso04_;
    float photonIso04_;

    bool mcMatched_;
    bool passConversionVeto_;
    bool isGsfCtfScPixChargeConsistent_;
    bool isPF_;
    
  };
}

#endif
