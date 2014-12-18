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

<<<<<<< HEAD
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
=======
    double chargedHadronIso(double dR=0.3) const {
>>>>>>> origin/70X_v2
      if( dR == 0.3){ return chargedHadronIso03_;}
      else if( dR == 0.4) { return chargedHadronIso04_;}
      else return -1.0;
    }
<<<<<<< HEAD
    float puChargedHadronIso(float dR=0.3) const {
=======
    double puChargedHadronIso(double dR=0.3) const {
>>>>>>> origin/70X_v2
      if( dR == 0.3) { return puChargedHadronIso03_;}
      else if( dR == 0.4) { return puChargedHadronIso04_;}
      else return -1.0;
    }
<<<<<<< HEAD
    float neutralHadronIso(float dR=0.3) const {
=======
    double neutralHadronIso(double dR=0.3) const {
>>>>>>> origin/70X_v2
      if( dR == 0.3) { return neutralHadronIso03_;}
      else if( dR == 0.4) { return neutralHadronIso04_;}
      else return -1.0;
    }
<<<<<<< HEAD
    float photonIso(float dR=0.3) const {
=======
    double photonIso(double dR=0.3) const {
>>>>>>> origin/70X_v2
      if( dR == 0.3) { return photonIso03_;}
      else if( dR == 0.4) { return photonIso04_;}
      else return -1.0;
    }

<<<<<<< HEAD
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

    float electronID(const std::string& name) const;
    float electronID(const char* name) const { return electronID( std::string(name) );}
    void setElectronIDs(const std::vector<pat::Electron::IdPair> & ids) { electronIDs_ = ids; }

  private:

    std::vector<pat::Electron::IdPair> electronIDs_;
    //std::vector<std::string> electronIDNames_;

    float relIso03_;
    float relIso04_;

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
=======
    double absIso(double dR=0.3, double dBetaFactor=0.5) const{

      if(dBetaFactor>0 && puChargedHadronIso(dR)<0) return -1;

      double neutralIso = neutralHadronIso(dR) + photonIso(dR);
      double corNeutralIso = neutralIso - dBetaFactor * puChargedHadronIso(dR);

      double charged = chargedHadronIso(dR);

      return charged + ( corNeutralIso>0 ? corNeutralIso : 0 ) ;
    }

    double relIso(double dR=0.3, double dBetaFactor=0.5) const{
      double abs = absIso(dR, dBetaFactor)/this->pt();
      return abs >=0 ? abs : -1;
    }


    bool mcMatched() const { return mcMatched_; }

    float mva() const { return mva_; }

    void setChargedHadronIso03(double i) { chargedHadronIso03_ = i; }
    void setPUChargedHadronIso03(double i) { puChargedHadronIso03_ = i; }
    void setNeutralHadronIso03(double i) { neutralHadronIso03_ = i; }
    void setPhotonIso03(double i) { photonIso03_ = i; }

    void setChargedHadronIso04(double i) { chargedHadronIso04_ = i; }
    void setPUChargedHadronIso04(double i) { puChargedHadronIso04_ = i; }
    void setNeutralHadronIso04(double i) { neutralHadronIso04_ = i; }
    void setPhotonIso04(double i) { photonIso04_ = i; }

    void setMCMatched(bool m) { mcMatched_ = m; }
    void setMva(float m) { mva_ = m; } 
  
  private:

    double chargedHadronIso03_;
    double puChargedHadronIso03_;
    double neutralHadronIso03_;
    double photonIso03_;

    double chargedHadronIso04_;
    double puChargedHadronIso04_;
    double neutralHadronIso04_;
    double photonIso04_;

    bool mcMatched_;
    float mva_;
>>>>>>> origin/70X_v2

    bool mcMatched_;
    bool passConversionVeto_;
    bool isGsfCtfScPixChargeConsistent_;
    bool isPF_;
    
  };
}

#endif
