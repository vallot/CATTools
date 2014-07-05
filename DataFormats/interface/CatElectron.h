#ifndef CatElectron_h
#define CatElectron_h

#include <map>
#include <iostream>
#include <TLorentzVector.h>


#include "../interface/CatLepton.h"

using namespace std;

namespace cat
{
	class CatElectron : public CatLepton
	{
    
	public:
		CatElectron() :
		    CatLepton(),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(9999.),
		    fBrem_(-9999.),
		    nBrems_(-9999),
		    Dist_(9999.),
		    DCot_(9999.),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;}
                CatElectron(const CatLepton& l) :
		    CatLepton(l),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(9999.),
		    fBrem_(-9999.),
		    nBrems_(-9999),
		    Dist_(9999.),
		    DCot_(9999.),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;} 
		CatElectron(const CatElectron& e) :
		    CatLepton(e),
		    trackerDrivenSeed_(e.trackerDrivenSeed_),
		    ecalDrivenSeed_(e.ecalDrivenSeed_),
		    ecalDrivenMomentum_(e.ecalDrivenMomentum_),
		    eSuperClusterOverPin_(e.eSuperClusterOverPin_),
		    eEleClusterOverPout_(e.eEleClusterOverPout_),
		    eSeedClusterOverPout_(e.eSeedClusterOverPout_),
		    deltaEtaIn_(e.deltaEtaIn_),
		    deltaEtaOut_(e.deltaEtaOut_),
		    deltaPhiIn_(e.deltaPhiIn_),
		    deltaPhiOut_(e.deltaPhiOut_),
		    deltaPhiSuperClusterTrackAtCalo_(e.deltaPhiSuperClusterTrackAtCalo_),
		    deltaEtaSuperClusterTrackAtCalo_(e.deltaEtaSuperClusterTrackAtCalo_),
		    ioEmIoP_(e.ioEmIoP_),
		    ioEmIoPgsf_(e.ioEmIoPgsf_),
		    pixelLayersWithMeasurement_(e.pixelLayersWithMeasurement_),
		    stripLayersWithMeasurement_(e.stripLayersWithMeasurement_),
		    nValidHits_(e.nValidHits_),
		    missingHits_(e.missingHits_),
		    normalizedChi2_(e.normalizedChi2_),
		    normalizedChi2gsf_(e.normalizedChi2gsf_),
		    superClusterRawEnergy_(e.superClusterRawEnergy_),
		    superClusterEta_(e.superClusterEta_),
		    preshowerEnergy_(e.preshowerEnergy_),
		    sigmaIetaIeta_(e.sigmaIetaIeta_),
		    sigmaIphiIphi_(e.sigmaIphiIphi_),
		    sigmaIetaIphi_(e.sigmaIetaIphi_),
		    e1x5_(e.e1x5_),
		    e5x5_(e.e5x5_),
		    hcalDepth1OverEcal_(e.hcalDepth1OverEcal_),
		    hcalDepth2OverEcal_(e.hcalDepth2OverEcal_),
		    etaWidth_(e.etaWidth_),
		    phiWidth_(e.phiWidth_),
		    r9_(e.r9_),
		    fBrem_(e.fBrem_),
		    nBrems_(e.nBrems_),
		    Dist_(e.Dist_),
		    DCot_(e.DCot_),
		    passConversion_(e.passConversion_),
		    mvaTrigId_(e.mvaTrigId_),
		    mvaNonTrigId_(e.mvaNonTrigId_)
		    {;}
    
		CatElectron(Double_t px, Double_t py, Double_t pz, Double_t e) :
		    CatLepton(px,py,pz,e),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(-9999.),
		    fBrem_(-9999.),
		    nBrems_(-9999),
		    Dist_(9999),
		    DCot_(9999),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;}
    
		CatElectron(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
		    CatLepton(px,py,pz,e,vtx_x,vtx_y,vtx_z),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(-9999.),
		    fBrem_(-9999.),
		    nBrems_(-9999),
		    Dist_(9999),
		    DCot_(9999),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;}
    
		CatElectron(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Int_t charge) :
		    CatLepton(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(-9999.),
		    fBrem_(-9999.),
		    nBrems_(-9999),
		    Dist_(9999),
		    DCot_(9999),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;}
    
		CatElectron(const TLorentzVector &momentum) :
		    CatLepton(momentum),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(-9999.),
		    fBrem_(-9999.),
		    nBrems_(9999),
		    Dist_(9999),
		    DCot_(9999),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;}
	    
                CatElectron(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
		    CatLepton(momentum, vertex, type, charge),
		    trackerDrivenSeed_(false),
		    ecalDrivenSeed_(false),
		    ecalDrivenMomentum_(),
		    eSuperClusterOverPin_(-9999.),
		    eEleClusterOverPout_(-9999.),
		    eSeedClusterOverPout_(-9999.),
		    deltaEtaIn_(-9999.),
		    deltaEtaOut_(-9999.),
		    deltaPhiIn_(-9999.),
		    deltaPhiOut_(-9999.),
		    deltaPhiSuperClusterTrackAtCalo_(-9999.),
		    deltaEtaSuperClusterTrackAtCalo_(-9999.),
		    ioEmIoP_(-9999.),
		    ioEmIoPgsf_(-9999.),
		    pixelLayersWithMeasurement_(-9999),
		    stripLayersWithMeasurement_(-9999),
		    nValidHits_(-9999),
		    missingHits_(-9999),
		    normalizedChi2_(9999.),
		    normalizedChi2gsf_(9999.),
		    superClusterRawEnergy_(-9999.),
		    superClusterEta_(-9999.),
		    preshowerEnergy_(-9999.),
		    sigmaIetaIeta_(-9999.),
		    sigmaIphiIphi_(-9999.),
		    sigmaIetaIphi_(-9999.),
		    e1x5_(-9999.),
		    e5x5_(-9999.),
		    hcalDepth1OverEcal_(-9999.),
		    hcalDepth2OverEcal_(-9999.),
		    etaWidth_(-9999.),
		    phiWidth_(-9999.),
		    r9_(-9999.),
		    fBrem_(-9999.),
		    nBrems_(9999),
		    Dist_(9999),
		    DCot_(9999),
		    passConversion_(false),
		    mvaTrigId_(-9999.),
		    mvaNonTrigId_(-9999.)
		    {;}
		    
		~CatElectron() {;}
	    
	    
	public:
	        virtual TString typeName() const { return "CatElectron"; }
	  
	        Bool_t isEcalDrivenSeed() const { return ecalDrivenSeed_; }
	        Bool_t isTrackerDrivenSeed() const { return trackerDrivenSeed_; }
	    
                TLorentzVector ecalDrivenMomentum() const {return ecalDrivenMomentum_; }
    
		Float_t eScOverP() const { return eSuperClusterOverPin_; }
                Float_t eEleClusterOverPout() const { return eEleClusterOverPout_; }
                Float_t eSeedClusterOverPout() const { return eSeedClusterOverPout_; }
		Float_t deltaEtaIn() const { return deltaEtaIn_; }
		Float_t deltaEtaOut() const { return deltaEtaOut_; }
		Float_t deltaPhiIn() const { return deltaPhiIn_; }
		Float_t deltaPhiOut() const { return deltaPhiOut_; }
		Float_t deltaPhiScTrkOut() const { return deltaPhiSuperClusterTrackAtCalo_; }
		Float_t deltaEtaScTrkOut() const { return deltaEtaSuperClusterTrackAtCalo_; }
                Float_t ioEmIoP() const { return ioEmIoP_; }
                Float_t ioEmIoPgsf() const { return ioEmIoPgsf_; }
    
		Int_t trackPixelLayersWithMeasurement() const { return pixelLayersWithMeasurement_; }
		Int_t trackStripLayersWithMeasurement() const { return stripLayersWithMeasurement_; }
                Int_t trackerLayersWithMeasurement() const { return pixelLayersWithMeasurement_+stripLayersWithMeasurement_; }
                Int_t trackNValidHits() const { return nValidHits_; }
		Int_t missingHits() const { return missingHits_; }
		Float_t trackNormalizedChi2() const { return normalizedChi2_; }
                Float_t gsfTrackNormalizedChi2() const { return normalizedChi2gsf_; }
    
		Float_t superClusterRawEnergy() const { return superClusterRawEnergy_; }
		Float_t superClusterEta() const { return superClusterEta_; }
		Float_t preshowerEnergy() const { return preshowerEnergy_; }
		Float_t sigmaIEtaIEta() const { return sigmaIetaIeta_; }
                Float_t sigmaIPhiIPhi() const { return sigmaIphiIphi_; }
                Float_t sigmaIEtaIPhi() const { return sigmaIetaIphi_; }
		Float_t e1x5() const { return e1x5_; }
		Float_t e5x5() const { return e5x5_; }
		Float_t hadronicOverEm() const { return (hcalDepth1OverEcal_ + hcalDepth2OverEcal_); }
		Float_t hadronicDepth1OverEm() const { return hcalDepth1OverEcal_; }
		Float_t hadronicDepth2OverEm() const { return hcalDepth2OverEcal_; }
                Float_t etaWidth() const { return etaWidth_; }
                Float_t phiWidth() const { return phiWidth_; }
                Float_t r9() const { return r9_; }
    
		Float_t fbrem() const { return fBrem_; }
		Int_t numberOfBrems() const { return nBrems_; }
		Float_t Dist() const { return Dist_; }
		Float_t DCot() const { return DCot_; }
                Bool_t passConversion() const { return passConversion_; }
                Float_t mvaTrigId() const { return mvaTrigId_; }
                Float_t mvaNonTrigId() const { return mvaNonTrigId_; }
    
		//setters
		void setEcalSeeding(Bool_t isEcal){ ecalDrivenSeed_ = isEcal; }
                void setEcalDrivenMomentum(TLorentzVector ecalDrivenMomentum){ ecalDrivenMomentum_ = ecalDrivenMomentum; }
		void setTrackerSeeding(Bool_t isTracker){ trackerDrivenSeed_ = isTracker; }
		void setDeltaEtaIn(Float_t deltaEtaIn) { deltaEtaIn_ = deltaEtaIn; }
		void setDeltaEtaOut(Float_t deltaEtaOut) { deltaEtaOut_ = deltaEtaOut; }
		void setDeltaEtaSuperClusterTrackAtCalo(Float_t x) { deltaEtaSuperClusterTrackAtCalo_ = x; }
		void setDeltaPhiIn(Float_t deltaPhiIn) { deltaPhiIn_ = deltaPhiIn; }
		void setDeltaPhiOut(Float_t deltaPhiOut) { deltaPhiOut_ = deltaPhiOut; }
		void setDeltaPhiSuperClusterTrackAtCalo(Float_t x) { deltaPhiSuperClusterTrackAtCalo_ = x; }
		void setEnergySuperClusterOverP(Float_t x) { eSuperClusterOverPin_ = x; }
                void setEnergyEleClusterOverPout(Float_t x) { eEleClusterOverPout_ = x; }
                void setEnergySeedClusterOverPout(Float_t x) { eSeedClusterOverPout_ = x; }
                void setIoEmIoP(Float_t x) { ioEmIoP_ = x; }
                void setIoEmIoPgsf(Float_t x) { ioEmIoPgsf_ = x; }
    
		void setTrackMissingHits(Int_t x) { missingHits_ = x; }
		void setTrackNormalizedChi2(Float_t x) { normalizedChi2_ = x; }
                void setGsfTrackNormalizedChi2(Float_t x) { normalizedChi2gsf_ = x; }
		void setPixelLayersWithMeasurement(Int_t x) { pixelLayersWithMeasurement_ = x; }
		void setStripLayersWithMeasurement(Int_t stripLayersWithMeasurement) { stripLayersWithMeasurement_ = stripLayersWithMeasurement; }
                void setNValidHits(Int_t nHits) { nValidHits_ = nHits; }
    
		void setPreshowerEnergy(Float_t x) { preshowerEnergy_ = x; }
		void setSuperClusterRawEnergy(Float_t x) { superClusterRawEnergy_ = x; }
		void setSuperClusterEta(Float_t x) { superClusterEta_ = x; }
    
		void setE1x5(Float_t e1x5) { e1x5_ = e1x5; }
		void setE5x5(Float_t e5x5) { e5x5_ = e5x5; }
		void setHoverEDepth1(Float_t HoE1) { hcalDepth1OverEcal_ = HoE1; }
		void setHoverEDepth2(Float_t HoE2) { hcalDepth2OverEcal_ = HoE2; }
		void setSigmaIetaIeta(Float_t sieie) { sigmaIetaIeta_ = sieie; }
                void setSigmaIphiIphi(Float_t sipip) { sigmaIphiIphi_ = sipip; }
                void setSigmaIetaIphi(Float_t sieip) { sigmaIetaIphi_ = sieip; }
                void setEtaWidth(Float_t etaWidth) { etaWidth_ = etaWidth; }
                void setPhiWidth(Float_t phiWidth) { phiWidth_ = phiWidth; }
                void setR9(Float_t r9) { r9_ = r9; }
    
		void setFbrem(Float_t f) { fBrem_ = f; }
                void setNBrems(Int_t n) { nBrems_ = n; }
		void setDist(Float_t dist) { Dist_ = dist; }
		void setDCot(Float_t dcot) { DCot_ = dcot; }
                void setPassConversion(Bool_t pass) { passConversion_ = pass; }
                void setMvaTrigId(Float_t id) { mvaTrigId_ = id; }
                void setMvaNonTrigId(Float_t id) { mvaNonTrigId_ = id; }
		
		friend std::ostream& operator<< (std::ostream& stream, const CatElectron& electron) {
			stream << "CatElectron - Charge=" << electron.charge() << " (Et,eta,phi)=("<< electron.Et() <<","<< electron.Eta() <<","<< electron.Phi() << ")"
      << " vertex(x,y,z)=("<< electron.vx() <<","<< electron.vy() <<","<< electron.vz() << ")";
			return stream;
      };
      
      
    private:
      Bool_t trackerDrivenSeed_;
      Bool_t ecalDrivenSeed_;
      TLorentzVector ecalDrivenMomentum_;        // ecal driven momentum, equivalent to gsf electron momentum.
      Float_t eSuperClusterOverPin_;             // the supercluster energy / track momentum at the PCA to the beam spot
      Float_t eEleClusterOverPout_;              // the electron cluster energy / track momentum at calo extrapolated from the outermost track state
      Float_t eSeedClusterOverPout_;             // the seed cluster energy / track momentum at calo extrapolated from the outermost track state
      Float_t deltaEtaIn_;                       // the supercluster eta - track eta position at calo extrapolated from innermost track state
      Float_t deltaEtaOut_;                      // the seed cluster eta - track eta position at calo extrapolated from the outermost track state
      Float_t deltaPhiIn_;                       // the supercluster phi - track phi position at calo extrapolated from the innermost track state
      Float_t deltaPhiOut_;                      // the seed cluster phi - track phi position at calo extrapolated from the outermost track state
      Float_t deltaPhiSuperClusterTrackAtCalo_;  // the electron cluster phi - track phi position at calo extrapolated from the outermost track state
      Float_t deltaEtaSuperClusterTrackAtCalo_;  // the electron cluster eta - t
      Float_t ioEmIoP_;                          // (1.0/(ele.superCluster()->energy())) - (1.0 / ele.p())
      Float_t ioEmIoPgsf_;                       // (1.0/(ele.superCluster()->energy())) - (1.0 / ele.gsfTrack()->p())
      
      //TrackProperties=====================================
      Int_t pixelLayersWithMeasurement_;         // number of pixel layers with at least one valid hit
      Int_t stripLayersWithMeasurement_;         // number of strip layers with at least one valid hit
      Int_t nValidHits_;                         // number of valid hits
      
      // In the standard PAT configuration, dB and edB are calculated wrt the primary vertex
      // If this was not the case, dB is calculated wrt the beamspot and edb = -1 all the time
      //Float_t dB_;                             // dB from PAT muon
      //Float_t dBError_;                        // dBError from PAT muon
      Int_t missingHits_;                        // Conversion Rejection: number of missing hits near beginning of track (also rejects really bad tracks)
      Float_t normalizedChi2_;                   // chi-squared divided by n.d.o.f. of track fit
      Float_t normalizedChi2gsf_;                // chi2 / ndf from gsfTrack
      
      //SuperClusterProperties ===============================
      Float_t superClusterRawEnergy_;
      Float_t superClusterEta_;
      Float_t preshowerEnergy_;
      
      //ShowerShape===========================================
      Float_t sigmaIetaIeta_;                    // weighted cluster rms along eta and inside 5x5 (new, Xtal eta)
      Float_t sigmaIphiIphi_;
      Float_t sigmaIetaIphi_;
      Float_t e1x5_;                             // energy inside 1x5 in etaxphi around the seed Xtal
      Float_t e5x5_;                             // energy inside 5x5 in etaxphi around the seed Xtal
      Float_t hcalDepth1OverEcal_ ;              // hcal over ecal seed cluster energy using first hcal depth (hcal is energy of towers within dR=015)
      Float_t hcalDepth2OverEcal_ ;              // hcal over ecal seed cluster energy using 2nd hcal depth (hcal is energy of towers within dR=015)
      Float_t etaWidth_;
      Float_t phiWidth_;
      Float_t r9_;
      
      // Electron classification && fBrem ====================
      Float_t fBrem_;                            // brem fraction from gsf fit: (track momentum in - track momentum out) / track momentum in
      Int_t   nBrems_;                           // number of basic clusters inside the supercluster - 1
      Float_t Dist_;                             // distance to the conversion partner
      Float_t DCot_;                             // difference of cot(angle) with the conversion partner track
      Bool_t  passConversion_;                   // boolean to flag converted candidates
      Float_t mvaTrigId_;                        // MVA value
      Float_t mvaNonTrigId_;                     // MVA value
      
      ClassDef (CatElectron,13);
      };
}

#endif
