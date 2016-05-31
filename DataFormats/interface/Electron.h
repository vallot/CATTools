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

    
    float shiftedEn() const { if (this->isEB()) return 0.006; else return 0.015; }
    float shiftedEnDown() const {return 1-shiftedEn();}
    float shiftedEnUp() const {return  1+shiftedEn();}
    
    int snuID() const {return snuID_;}
    bool isTrigMVAValid() const { return isTrigMVAValid_; }
    
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

    void setSNUID(int id) {snuID_ = id;}
    void setTrigMVAValid(bool val) { isTrigMVAValid_ = val; }
    void setIpSignficance(float ipsig) {ipsig_ = ipsig;}


    float scaleFactor(const std::string& name, int sign = 0) const;
    
  private:

    std::vector<pat::Electron::IdPair> electronIDs_;

    float relIso03_;
    float relIso04_;
    float ipsig_;
    float scEta_;
    bool passConversionVeto_;
    bool isGsfCtfScPixChargeConsistent_;
    bool isEB_;
    
    int snuID_;
    bool isTrigMVAValid_;
  };
}

#endif
