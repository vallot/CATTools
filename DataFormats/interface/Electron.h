#ifndef CATTools_Electron_H
#define CATTools_Electron_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

// Define typedefs for convenience
namespace cat {
  class Electron;
  typedef std::vector<Electron>              ElectronCollection;
  typedef edm::Ref<ElectronCollection>       ElectronRef;
  typedef edm::RefVector<ElectronCollection> ElectronRefVector;
}

namespace cat {

  class Electron : public Lepton{
  public:
    Electron();
    Electron(const reco::LeafCandidate & aElectron);
    virtual ~Electron();

    float electronID(const std::string& name) const;
    float electronID(const char* name) const { return electronID( std::string(name) );}
    bool isVeto() const {return electronID("veto");}
    bool isMediumMVA() const {return electronID("wp90");}
    bool isTightMVA() const {return electronID("wp80");}

    float relIso(float dR=0.3) const {
      if( dR < 0.35) return relIso03_;
      else return relIso04_;
    }
    
    float scEta() const { return scEta_; }
    bool passConversionVeto() const { return passConversionVeto_; }
    bool isGsfCtfScPixChargeConsistent() const{ return isGsfCtfScPixChargeConsistent_; }
    bool isEB() const{ return isEB_; }

    float ipsignificance() const { return ipsig_;}

    float smearedScale() const { return smearedScale_[0]; }
    float shiftedEn() const { return 0.0; } //Put dummy number
    float shiftedEnDown( int i = 0 ) const {
      if ( i == 0 ) return smearedScale_[2];   // userFloat("energyScaleDown")
      else          return smearedScale_[4]; } // userFloat("energySigmaDown")
    float shiftedEnUp( int i = 0 ) const {
      if ( i == 0 ) return smearedScale_[1];   // userFloat("energyScaleUp")
      else          return smearedScale_[3]; } // userFloat("energySigmaUp")
    
    void setElectronIDs(const std::vector<pat::Electron::IdPair> & ids) { electronIDs_ = ids; }
    void setElectronID(pat::Electron::IdPair ids) { electronIDs_.push_back(ids); }

    void setrelIso(double dR, double chIso, double nhIso, double phIso, double AEff, double rhoIso, double ecalpt)
    {
      float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;
      if( dR < 0.35) relIso03_ = relIso;
      else  relIso04_ = relIso;
    }
    
    void setscEta(float i) { scEta_ = i; }
    void setPassConversionVeto(bool i) {  passConversionVeto_ = i; }
    void setIsGsfCtfScPixChargeConsistent(bool d) { isGsfCtfScPixChargeConsistent_ = d ; }
    void setIsEB(bool d) { isEB_ = d ; }

    void setIpSignficance(float ipsig) {ipsig_ = ipsig;}

    void setSmearedScale(const float scale) {
      smearedScale_[0] = scale; } //userFloat("ecalTrkEnergyPostCorr") / energy()
    void setSmearedScaleUnc(const float up1, const float dn1, const float up2, const float dn2) {
      smearedScale_[1] = up1;   // userFloat("energyScaleUp")
      smearedScale_[2] = dn1;   // userFloat("energyScaleDown")
      smearedScale_[3] = up2;   // userFloat("energySigmaUp")
      smearedScale_[4] = dn2; } // userFloat("energySigmaDown")

  private:

    std::vector<pat::Electron::IdPair> electronIDs_;

    float smearedScale_[5];
    float relIso03_;
    float relIso04_;
    float ipsig_;
    float scEta_;
    bool passConversionVeto_;
    bool isGsfCtfScPixChargeConsistent_;
    bool isEB_;
    
  };
}

#endif
