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

    int GlobalCharge() const { return GlobalCharge_; }
		float GlobalEta() const { return GlobalEta_; }
		float GlobalPt() const { return GlobalPt_; }
		float GlobalPhi() const { return GlobalPhi_; }
  	float Energy() const { return Energy_; }
  
    float PtError() const { return PtError_; }
    float EtaError() const { return EtaError_; }
    
    float PFIsoR03ChargedHadron() const { return PFIsoR03ChargedHadron_; }
    float PFIsoR03NeutralHadron() const { return PFIsoR03NeutralHadron_; }
    float PFIsoR03Photon() const { return PFIsoR03Photon_; }
    float PFIsoR04ChargedHadron() const { return PFIsoR04ChargedHadron_; }
    float PFIsoR04NeutralHadron() const { return PFIsoR04NeutralHadron_; }
    float PFIsoR04Photon() const { return PFIsoR04Photon_; }
    float EcalVetoIso() const { return EcalVetoIso_; }
    float HcalVetoIso() const { return HcalVetoIso_; }
    float PFIsoR03PU() const { return PFIsoR03PU_; }
    float PFIsoR04PU() const { return PFIsoR04PU_; }
    
		float TrkVx() const { return TrkVx_; }
		float TrkVy() const { return TrkVy_; }
		float TrkVz() const { return TrkVz_; }

		float MatchedGenParticlePt() const { return MatchedGenParticlePt_; }
		float MatchedGenParticleEta() const { return MatchedGenParticleEta_; }
		float MatchedGenParticlePhi() const { return MatchedGenParticlePhi_; }
		float TrkD0() const { return TrkD0_; }
		float TrkD0Error() const { return TrkD0Error_; }	

		float GlobalChi2() const { return GlobalChi2_; }
		float PrimaryVertexDXY() const { return PrimaryVertexDXY_; }
		float PrimaryVertexDXYError() const { return PrimaryVertexDXYError_; }
		
		int GlobalTrkValidHits() const { return GlobalTrkValidHits_; }
		int TrkPixelHits() const { return TrkPixelHits_; }
		int StationMatches() const { return StationMatches_; }
		int TrackLayersWithMeasurement() const { return TrackLayersWithMeasurement_; }
		int IsPF() const { return IsPF_; }
		int IsGlobal() const { return IsGlobal_; }
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

		int TrackerCharge() const { return TrackerCharge_; }
		
		float VtxDistXY() const { return VtxDistXY_; }
		int BestTrackVtxIndex() const { return BestTrackVtxIndex_; }
		float BestTrackVtxDistZ() const { return BestTrackVtxDistZ_; }
		float BestTrackVtxDistXY() const { return BestTrackVtxDistXY_; }


		///////////////////
		// set functions //
		///////////////////
    
    void setGlobalCharge(int i) { GlobalCharge_ = i; }
		void setGlobalEta(int f) { GlobalEta_ = f; }
		void setGlobalPt(float f) { GlobalPt_ = f;}
		void setGlobalPhi(float f) { GlobalPhi_ = f;}
   	void setEnergy(float f) { Energy_ = f; }
 
    void setPtError(float f) { PtError_ = f; }
    void setEtaError(float f) { EtaError_ = f; }
    
    void setPFIsoR03ChargedHadron(float f) { PFIsoR03ChargedHadron_ = f; }
    void setPFIsoR03NeutralHadron(float f) { PFIsoR03NeutralHadron_ = f; }
    void setPFIsoR03Photon(float f) { PFIsoR03Photon_ = f; }
    void setPFIsoR04ChargedHadron(float f) { PFIsoR04ChargedHadron_ = f; }
    void setPFIsoR04NeutralHadron(float f) { PFIsoR04NeutralHadron_ = f; }
    void setPFIsoR04Photon(float f) { PFIsoR04Photon_ = f; }
    void setEcalVetoIso(float f) { EcalVetoIso_ = f; }
    void setHcalVetoIso(float f) { HcalVetoIso_ = f; }
    void setPFIsoR03PU(float f) { PFIsoR03PU_ = f; }
    void setPFIsoR04PU(float f) { PFIsoR04PU_ = f; }
    
		void setTrkVx(float f) { TrkVx_ = f; }
		void setTrkVy(float f) { TrkVy_ = f; }
		void setTrkVz(float f) { TrkVz_ = f; }

		void setMatchedGenParticlePt(float f) { MatchedGenParticlePt_ = f; }		
		void setMatchedGenParticleEta(float f) { MatchedGenParticleEta_ = f; }
		void setMatchedGenParticlePhi(float f) { MatchedGenParticlePhi_ = f; }
		void setTrkD0(float f) { TrkD0_ = f; }
		void setTrkD0Error(float f) { TrkD0Error_ = f; }

		void setGlobalChi2(float f) { GlobalChi2_ = f; }
		void setPrimaryVertexDXY(float f) { PrimaryVertexDXY_ = f; }
		void setPrimaryVertexDXYError(float f) { PrimaryVertexDXYError_ = f; }

    void setGlobalTrkValidHits(int i) { GlobalTrkValidHits_ = i; }
    void setTrkPixelHits(int i) { TrkPixelHits_ = i; }
    void setStationMatches(int i) { StationMatches_ = i; }
    void setTrackLayersWithMeasurement(int i) { TrackLayersWithMeasurement_ = i; }
    void setIsPF(int i) { IsPF_ = i; }
    void setIsGlobal(int i) { IsGlobal_ = i; }
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

		void setTrackerCharge(int i) { TrackerCharge_ = i; }

    void setVtxDistXY(float f) { VtxDistXY_ = f; }
    void setBestTrackVtxIndex(int i) { BestTrackVtxIndex_ = i; }
    void setBestTrackVtxDistZ(float f) { BestTrackVtxDistZ_ = f; }
    void setBestTrackVtxDistXY(float f) { BestTrackVtxDistXY_ = f; }

  private:

    int GlobalCharge_;
		float GlobalEta_;
		float GlobalPt_;
		float GlobalPhi_;
		float Energy_;   
 
    float PtError_;
    float EtaError_;
    
    float PFIsoR03ChargedHadron_;
    float PFIsoR03NeutralHadron_;
    float PFIsoR03Photon_;
    float PFIsoR04ChargedHadron_;
    float PFIsoR04NeutralHadron_;
    float PFIsoR04Photon_;
    float EcalVetoIso_;
    float HcalVetoIso_;
    float PFIsoR03PU_;
    float PFIsoR04PU_;

		float TrkVx_;
		float TrkVy_;
		float TrkVz_;

		float MatchedGenParticlePt_;
		float MatchedGenParticleEta_;
		float MatchedGenParticlePhi_;
		float TrkD0_;
		float TrkD0Error_;
		
		float GlobalChi2_;
		float PrimaryVertexDXY_;
		float PrimaryVertexDXYError_;
	
		int GlobalTrkValidHits_;
		int TrkPixelHits_;
		int StationMatches_;
		int TrackLayersWithMeasurement_;
		int IsPF_;
		int IsGlobal_;
		int IsTracker_;

		float CocktailPt_, CocktailEta_, CocktailPhi_, CocktailGlobalChi2_, CocktailTrkVtxDXY_, CocktailTrkVtxDZ_;
		int CocktailCharge_;
		
		float MuonSpecPt_, MuonSpecEta_, MuonSpecPhi_, MuonSpecE_;
		int MuonSpecCharge_;

		int TrackerCharge_;

		float VtxDistXY_, BestTrackVtxIndex_, BestTrackVtxDistZ_, BestTrackVtxDistXY_;
  };
}

#endif
