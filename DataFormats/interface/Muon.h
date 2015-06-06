#ifndef CATTools_Muon_H
#define CATTools_Muon_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// Define typedefs for convenience
namespace cat {
  class Muon;
  typedef std::vector<Muon>              MuonCollection;
  typedef edm::Ref<MuonCollection>       MuonRef;
  typedef edm::RefVector<MuonCollection> MuonRefVector;
}

namespace cat {

  class Muon : public Particle{
  public:
    Muon();
    Muon(const reco::LeafCandidate & aMuon); 
    virtual ~Muon();

    /// cat default variables ///
    
    float relIso(float dR=0.3 ) const {
      if( dR < 0.35) return relIso03_;
      else return relIso04_;
    }
    bool isGlobalMuon() const { return isGlobalMuon_; }
    bool isPFMuon() const { return isPFMuon_; }
    bool isTightMuon() const { return isTightMuon_; }
    bool isLooseMuon() const { return isLooseMuon_; } 
    bool isSoftMuon() const { return isSoftMuon_; } 

    bool mcMatched() const { return mcMatched_; }

    float normalizedChi2() const { return normalizedChi2_; }
    int numberOfValidHits() const { return numberOfValidHits_; }
    int numberOfValidMuonHits() const { return numberOfValidMuonHits_; }
    int numberOfMatchedStations() const { return numberOfMatchedStations_; } 
    int numberOfValidPixelHits() const { return numberOfValidPixelHits_; }
    int trackerLayersWithMeasurement() const { return trackerLayersWithMeasurement_; }
 
    float dxy() const { return dxy_; }
    float dz() const { return dz_; }

    void setrelIso(float dR, double chIso, double nhIso, double phIso, double puIso, double pt) {
      float relIso = ( chIso + std::max( 0.0, nhIso + phIso - 0.5 * puIso) ) / pt;
      if( dR < 0.35) relIso03_ = relIso; 
      else  relIso04_ = relIso; 
    }
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

    void setChargedHadronIso03(float i) { chargedHadronIso03_ = i; }
    void setPUChargedHadronIso03(float i) { puChargedHadronIso03_ = i; }
    void setNeutralHadronIso03(float i) { neutralHadronIso03_ = i; }
    void setPhotonIso03(float i) { photonIso03_ = i; }

    void setChargedHadronIso04(float i) { chargedHadronIso04_ = i; }
    void setPUChargedHadronIso04(float i) { puChargedHadronIso04_ = i; }
    void setNeutralHadronIso04(float i) { neutralHadronIso04_ = i; }
    void setPhotonIso04(float i) { photonIso04_ = i; }

    void setIsGlobalMuon(bool d) { isGlobalMuon_ = d; }
    void setIsPFMuon(bool d) { isPFMuon_ = d; }
    void setIsTightMuon(bool d) { isTightMuon_ = d; }
    void setIsLooseMuon(bool d) { isLooseMuon_ = d; }
    void setIsSoftMuon(bool d) { isSoftMuon_ = d; }
     
    void setMCMatched(bool m) { mcMatched_ = m; }
    
    void setNormalizedChi2(float d) { normalizedChi2_ = d; }
    void setNumberOfValidHits(int i) { numberOfValidHits_ = i; }
    void setNumberOfValidMuonHits(int i) { numberOfValidMuonHits_ = i; }
    void setNumberOfMatchedStations(int i) { numberOfMatchedStations_ = i; }
    void setNumberOfValidPixelHits(int i) { numberOfValidPixelHits_ = i; }
    void setTackerLayersWithMeasurement(int i) { trackerLayersWithMeasurement_ = i; }
 
    void setDz(float d) { dz_ = d; }
    void setDxy(float d) { dxy_ = d; }

    void setShiftedEnDown(float f) { shiftedEnDown_ = f;}
    void setShiftedEnUp(float f) { shiftedEnUp_ = f;}
    float shiftedEnDown() {return  shiftedEnDown_;}
    float shiftedEnUp()   {return  shiftedEnUp_;}

    /// SKTree variables ///
    
    float EcalVetoIso() const { return EcalVetoIso_; }
    float HcalVetoIso() const { return HcalVetoIso_; }
    float TrkVx() const { return TrkVx_; }
    float TrkVy() const { return TrkVy_; }
    float TrkVz() const { return TrkVz_; }
    float MatchedGenParticlePt() const { return MatchedGenParticlePt_; }
    float MatchedGenParticleEta() const { return MatchedGenParticleEta_; }
    float MatchedGenParticlePhi() const { return MatchedGenParticlePhi_; }
    int IsTracker() const { return IsTracker_; }
    float CocktailPt() const { return CocktailPt_; }
    float CocktailEta() const { return CocktailEta_; }
    float CocktailPhi() const { return CocktailPhi_; }
    float CocktailGlobalChi2() const { return CocktailGlobalChi2_; }
    float CocktailTrkVtxDXY() const { return CocktailTrkVtxDXY_; }
    float CocktailTrkVtxDZ() const { return CocktailTrkVtxDZ_; }
    int CocktailCharge() const { return CocktailCharge_; }
    float MuonSpecPt() const { return MuonSpecPt_; }
    float MuonSpecEta() const { return MuonSpecEta_; }
    float MuonSpecPhi() const { return MuonSpecPhi_; }
    int MuonSpecCharge() const { return MuonSpecCharge_; }
    float MuonSpecE() const { return MuonSpecE_; }
    
    void setEcalVetoIso(float f) { EcalVetoIso_ = f; }
    void setHcalVetoIso(float f) { HcalVetoIso_ = f; }
    void setTrkVx(float f) { TrkVx_ = f; }
    void setTrkVy(float f) { TrkVy_ = f; }
    void setTrkVz(float f) { TrkVz_ = f; }
    void setMatchedGenParticlePt(float f) { MatchedGenParticlePt_ = f; }
    void setMatchedGenParticleEta(float f) { MatchedGenParticleEta_ = f; }
    void setMatchedGenParticlePhi(float f) { MatchedGenParticlePhi_ = f; }
    void setIsTracker(int i) { IsTracker_ = i; }
    void setCocktailPt(float f) {  CocktailPt_ = f; }
    void setCocktailEta(float f) {  CocktailEta_ = f; }
    void setCocktailPhi(float f) {  CocktailPhi_ = f; }
    void setCocktailGlobalChi2(float f) {  CocktailGlobalChi2_ = f; }
    void setCocktailTrkVtxDXY(float f) {  CocktailTrkVtxDXY_ = f; }
    void setCocktailTrkVtxDZ(float f) {  CocktailTrkVtxDZ_ = f; }
    void setCocktailCharge(int i) {  CocktailCharge_ = i; }
    void setMuonSpecPt(float f) {  MuonSpecPt_ = f; }
    void setMuonSpecEta(float f) {  MuonSpecEta_ = f; }
    void setMuonSpecPhi(float f) {  MuonSpecPhi_ = f; }
    void setMuonSpecCharge(int i) {  MuonSpecCharge_ = i; }
    void setMuonSpecE(float f) {  MuonSpecE_ = f; }

    
  private:

    /// cat default variables ///
    
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

    bool isGlobalMuon_; 
    bool isPFMuon_; 
    bool isTightMuon_; 
    bool isLooseMuon_; 
    bool isSoftMuon_;

    bool mcMatched_;

    float normalizedChi2_; 
    int numberOfValidHits_; 
    int numberOfValidMuonHits_; 
    int numberOfMatchedStations_; 
    int numberOfValidPixelHits_; 
    int trackerLayersWithMeasurement_; 

    float dz_;
    float dxy_;

    float shiftedEnDown_;
    float shiftedEnUp_;

    /// SKTree variables ///
    
    float EcalVetoIso_, HcalVetoIso_;
    float TrkVx, TrkVy_, TrkVz_;
    float MatchedGenParticlePt_, MatchedGenParticleEta_, MatchedGenParticlePhi_;
    int IsTracker_;
    float CocktailPt_, CocktailEta_, CocktailPhi_, CocktailGlobalChi2_, CocktailTrkVtxDXY_, CocktailTrkVtxDZ_;
    int CocktailCharge_;
    float MuonSpecPt_, MuonSpecEta_, MuonSpecPhi_, MuonSpecE_;
    int MuonSpecCharge_;
    
  };
}

#endif

