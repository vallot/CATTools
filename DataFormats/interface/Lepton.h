#ifndef CATTools_Lepton_H
#define CATTools_Lepton_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"

// Define typedefs for convenience
namespace cat {
  class Lepton;
  typedef std::vector<Lepton>              LeptonCollection;
  typedef edm::Ref<LeptonCollection>       LeptonRef;
  typedef edm::RefVector<LeptonCollection> LeptonRefVector;
}

namespace cat {

  class Lepton : public Particle{
  public:
    Lepton();
    Lepton(const reco::LeafCandidate & aLepton);
    virtual ~Lepton();

    bool isPF() const{ return isPF_; }
    bool isTight() const { return isTight_; }
    bool isMedium() const { return isMedium_; }
    bool isLoose() const { return isLoose_; }
    
    float dxy() const { return dxy_; }
    float dz() const { return dz_; }

    float chargedHadronIso(float dR=0.3) const {
      if( dR < 0.35) return chargedHadronIso03_;
      else return chargedHadronIso04_;
    }
    float puChargedHadronIso(float dR=0.3) const {
      if( dR < 0.35) return puChargedHadronIso03_;
      else return puChargedHadronIso04_;
    }
    float neutralHadronIso(float dR=0.3) const {
      if( dR < 0.35) return neutralHadronIso03_;
      else return neutralHadronIso04_;
    }
    float photonIso(float dR=0.3) const {
      if( dR < 0.35) return photonIso03_;
      else return photonIso04_;
    }
    float absIso(float dR=0.3, float dBetaFactor=0.5) const{
      if(dBetaFactor>0 && puChargedHadronIso(dR)<0) return -1;
      float neutralIso = neutralHadronIso(dR) + photonIso(dR);
      float corNeutralIso = neutralIso - dBetaFactor * puChargedHadronIso(dR);
      float charged = chargedHadronIso(dR);
      return charged + ( corNeutralIso>0 ? corNeutralIso : 0 ) ;
    }
    float relIso(float dR=0.3, float dBetaFactor = 0.5) const {
	float abs = absIso(dR, dBetaFactor)/this->pt();
	return abs >=0 ? abs : -1;
    }
    // to be undated with shifts on the fly!
    /* float shiftedEnDown() const {return  shiftedEnDown_;} */
    /* float shiftedEnUp() const {return  shiftedEnUp_;} */

    bool mcMatched() const { return mcMatched_; }

    void setIsPF(bool d) { isPF_ = d ; }
    void setIsTight(bool d) { isTight_ = d; }
    void setIsMedium(bool d) { isMedium_ = d; }
    void setIsLoose(bool d) { isLoose_ = d; }

    void setDz(float d) { dz_ = d; }
    void setDxy(float d) { dxy_ = d; }
    
    void setChargedHadronIso03(float i) { chargedHadronIso03_ = i; }
    void setPUChargedHadronIso03(float i) { puChargedHadronIso03_ = i; }
    void setNeutralHadronIso03(float i) { neutralHadronIso03_ = i; }
    void setPhotonIso03(float i) { photonIso03_ = i; }

    void setChargedHadronIso04(float i) { chargedHadronIso04_ = i; }
    void setPUChargedHadronIso04(float i) { puChargedHadronIso04_ = i; }
    void setNeutralHadronIso04(float i) { neutralHadronIso04_ = i; }
    void setPhotonIso04(float i) { photonIso04_ = i; }

    void setMCMatched(bool m) { mcMatched_ = m; }
    
  private:

    bool isPF_;
    bool isTight_;
    bool isMedium_;
    bool isLoose_;    

    float dz_;
    float dxy_;

    float chargedHadronIso03_;
    float puChargedHadronIso03_;
    float neutralHadronIso03_;
    float photonIso03_;

    float chargedHadronIso04_;
    float puChargedHadronIso04_;
    float neutralHadronIso04_;
    float photonIso04_;

    bool mcMatched_;
  };
}

#endif
