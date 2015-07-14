#ifndef CATTools_Electron_H
#define CATTools_Electron_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
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
    float scPhi() const { return scPhi_; }
    float scPt() const { return scPt_; }
    float scRawEnergy() const { return scRawEnergy_; }
    float dxy() const { return dxy_; }
    float dz() const { return dz_; }
    float rho() const { return rho_; }
    bool passConversionVeto() const { return passConversionVeto_; }
    bool isGsfCtfScPixChargeConsistent() const { return isGsfCtfScPixChargeConsistent_; }
    bool isGsfScPixChargeConsistent() const { return isGsfScPixChargeConsistent_; }
    bool isGsfCtfChargeConsistent() const { return isGsfCtfChargeConsistent_; }

    bool isEB() const { return isEB_; }
    bool isEE() const { return isEE_; }
    bool trackerDrivenSeed() const { return trackerDrivenSeed_; }
    bool ecalDrivenSeed() const { return ecalDrivenSeed_; }

    float dr03TkSumPt() const { return trkIsoDR03_; }
    float dr03EcalRecHitSumEt() const { return ecalIsoDR03_; }
    float dr03HcalTowerSumEt() const { return hcalIsoDR03_; }
    float dr04TkSumPt() const { return trkIsoDR04_; }
    float dr04EcalRecHitSumEt() const { return ecalIsoDR04_; }
    float dr04HcalTowerSumEt() const { return hcalIsoDR04_; }

    float deltaPhiSuperClusterTrackAtVtx() const { return deltaPhiTrkSC_; }
    float deltaEtaSuperClusterTrackAtVtx() const { return deltaEtaTrkSC_; }

    float sigmaIEtaIEta() const { return sigmaIEtaIEta_; }

    float hadronicOverEm() const { return hoe_; }
    float caloEnergy() const { return caloEnergy_; }
    float eSuperClusterOverP() const { return eSuperClusterOverP_; }
    
    int numberOfBrems() const { return numberOfBrems_; }
    float fbrem() const { return fbrem_; }

    float dB() const { return dB_; }
    float edB() const { return edB_; }
    
    bool isPF() const { return isPF_; }

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

    bool mcMatched() const { return mcMatched_; }

///////////////////////////////////

    void setrelIso(double dR, double chIso, double nhIso, double phIso, double AEff, double rhoIso, double ecalpt)
    {
      float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;
      if( dR < 0.35) relIso03_ = relIso; 
      else  relIso04_ = relIso; 
    }
    void setscEta(float i) { scEta_ = i; } 
    void setscPhi(float i) { scPhi_ = i; }
    void setscPt(float i) { scPt_ = i; }
    void setscRawEnergy(float i) { scRawEnergy_ = i; }
    void setdxy(float i) {  dxy_ = i; }
    void setdz(float i) {  dz_ = i; }
    void setPassConversionVeto(bool i) {  passConversionVeto_ = i; }
    void setIsGsfCtfScPixChargeConsistent(bool i) {  isGsfCtfScPixChargeConsistent_ = i; }
    void setIsGsfScPixChargeConsistent(bool i) {  isGsfScPixChargeConsistent_ = i; }
    void setIsGsfCtfChargeConsistent(bool i) {  isGsfCtfChargeConsistent_ = i; }

    void setIsEB(bool i) { isEB_ = i; }
    void setIsEE(bool i) { isEE_ = i; }
    void setTrackerDrivenSeed(bool i) { trackerDrivenSeed_ = i; }
    void setEcalDrivenSeed(bool i) { ecalDrivenSeed_ = i; }

    void setdr03TkSumPt(float i) { trkIsoDR03_ = i; }
    void setdr03EcalRecHitSumEt(float i) { ecalIsoDR03_ = i; }
    void setdr03HcalTowerSumEt(float i) { hcalIsoDR03_ = i; }
    void setdr04TkSumPt(float i) { trkIsoDR04_ = i; }
    void setdr04EcalRecHitSumEt(float i) { ecalIsoDR04_ = i; }
    void setdr04HcalTowerSumEt(float i) { hcalIsoDR04_ = i; }

    void setdeltaPhiSuperClusterTrackAtVtx(float i) { deltaPhiTrkSC_ = i; }
    void setdeltaEtaSuperClusterTrackAtVtx(float i) { deltaEtaTrkSC_ = i; }

    void setSigmaIEtaIEta(float i) { sigmaIEtaIEta_ = i; }
    void sethadronicOverEm(float i) { hoe_ = i; }
    void setCaloEnergy(float i) { caloEnergy_ = i; }
    void setESuperClusterOverP(float i) { eSuperClusterOverP_ = i; }

    void setNumberOfBrems(int i) { numberOfBrems_ = i; }
    void setFbrem(float i) { fbrem_ = i; }

    void setdB(float i) { dB_ = i; }
    void setedB(float i) { edB_ = i; }
    
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
    void setElectronID(pat::Electron::IdPair ids) { electronIDs_.push_back(ids); }

    void setShiftedEnDown(float f) { shiftedEnDown_ = f;}
    void setShiftedEnUp(float f) { shiftedEnUp_ = f;}
    float shiftedEnDown() const {return  shiftedEnDown_;}
    float shiftedEnUp() const {return  shiftedEnUp_;}

    //    const reco::GsfTrack & gsfTrack() const {return  gsfTrack_;}
    //    void setgsfTrack(const reco::GsfTrack & gsf) { gsfTrack_ = gsf;}    
      
  private:

    std::vector<pat::Electron::IdPair> electronIDs_;
    //    reco::GsfTrack gsfTrack_;
    
    bool isEB_;
    bool isEE_;
    bool trackerDrivenSeed_;
    bool ecalDrivenSeed_;

    float relIso03_;
    float relIso04_;

    float chargedHadronIso03_;
    float puChargedHadronIso03_;
    float neutralHadronIso03_;
    float photonIso03_;

    float chargedHadronIso04_;
    float puChargedHadronIso04_;
    float neutralHadronIso04_;
    float photonIso04_;

    float scEta_;
    float scPhi_;
    float scPt_;
    float scRawEnergy_;
    float dxy_;
    float dz_;
    float rho_;
    
    float trkIsoDR03_;
    float ecalIsoDR03_;
    float hcalIsoDR03_;
    float trkIsoDR04_;
    float ecalIsoDR04_;
    float hcalIsoDR04_;

    float deltaPhiTrkSC_;
    float deltaEtaTrkSC_;

    float sigmaIEtaIEta_;
    float hoe_;
    float caloEnergy_;
    float eSuperClusterOverP_;
    float trackVx_;
    float trackVy_;
    float trackVz_;

    int numberOfBrems_;
    float fbrem_;
    float dB_;
    float edB_;
    float trackPt_;
    float trackValidFractionOfHits_;

    int vtxIndex_;
    float vtxDistXY_;
    float vtxDistZ_;
    float leadVtxDistXY_;
    float leadVtxDistZ_;

    bool mcMatched_;
    bool passConversionVeto_;
    bool isGsfCtfScPixChargeConsistent_;
    bool isGsfScPixChargeConsistent_;
    bool isGsfCtfChargeConsistent_;
    bool isPF_;

    float shiftedEnDown_;
    float shiftedEnUp_;

  };
}

#endif
