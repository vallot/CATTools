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

    double chargedHadronIso(double dR=0.3) const {
      if( dR == 0.3){ return chargedHadronIso03_;}
      else if( dR == 0.4) { return chargedHadronIso04_;}
      else return -1.0;
    }
    double puChargedHadronIso(double dR=0.3) const {
      if( dR == 0.3) { return puChargedHadronIso03_;}
      else if( dR == 0.4) { return puChargedHadronIso04_;}
      else return -1.0;
    }
    double neutralHadronIso(double dR=0.3) const {
      if( dR == 0.3) { return neutralHadronIso03_;}
      else if( dR == 0.4) { return neutralHadronIso04_;}
      else return -1.0;
    }
    double photonIso(double dR=0.3) const {
      if( dR == 0.3) { return photonIso03_;}
      else if( dR == 0.4) { return photonIso04_;}
      else return -1.0;
    }

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

  };
}

#endif
