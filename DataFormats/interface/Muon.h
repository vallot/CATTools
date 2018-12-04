#ifndef CATTools_Muon_H
#define CATTools_Muon_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Define typedefs for convenience
namespace cat {
  class Muon;
  typedef std::vector<Muon>              MuonCollection;
  typedef edm::Ref<MuonCollection>       MuonRef;
  typedef edm::RefVector<MuonCollection> MuonRefVector;
}

namespace cat {

  class Muon : public Lepton{
  public:
    Muon();
    Muon(const reco::LeafCandidate & aMuon);
    virtual ~Muon();
    
    bool isGlobalMuon() const { return isGlobalMuon_; }
    bool isSoftMuon() const { return isSoftMuon_; }
    bool isPFMuon() const { return this->isPF(); }
    bool isTrackerMuon() const { return isTrackerMuon_; }
    bool isTightMuon() const { return this->isTight(); }
    bool isMediumMuon() const { return this->isMedium(); }
    bool isLooseMuon() const { return this->isLoose(); }

    float normalizedChi2() const { return normalizedChi2_; }
    int numberOfValidHits() const { return numberOfValidHits_; }
    int numberOfValidMuonHits() const { return numberOfValidMuonHits_; }
    int numberOfMatchedStations() const { return numberOfMatchedStations_; }
    int numberOfValidPixelHits() const { return numberOfValidPixelHits_; }
    int trackerLayersWithMeasurement() const { return trackerLayersWithMeasurement_; }

    float ipsignificance() const { return ipsig_;}


    float shiftedEn() const { if (this->pt() < 100) return 0.002; else return 0.05; }
    float shiftedEnDown() const {return 1-shiftedEn();}
    float shiftedEnUp() const {return  1+shiftedEn();}

    void setIsGlobalMuon(bool d) { isGlobalMuon_ = d; }
    void setIsSoftMuon(bool d) { isSoftMuon_ = d; }
    void setIsTrackerMuon(bool d) { isTrackerMuon_ = d; }

    void setNormalizedChi2(float d) { normalizedChi2_ = d; }
    void setNumberOfValidHits(int i) { numberOfValidHits_ = i; }
    void setNumberOfValidMuonHits(int i) { numberOfValidMuonHits_ = i; }
    void setNumberOfMatchedStations(int i) { numberOfMatchedStations_ = i; }
    void setNumberOfValidPixelHits(int i) { numberOfValidPixelHits_ = i; }
    void setTackerLayersWithMeasurement(int i) { trackerLayersWithMeasurement_ = i; }

    void setIpSignficance(float ipsig) {ipsig_ = ipsig;}

    float scaleFactor(const std::string& name, int sign = 0) const;

    /// DataFormats/MuonReco/interface/Muon.h
    /// ====================== STANDARD SELECTORS ===========================
    ///
    enum Selector {
      CutBasedIdLoose        = 1UL<< 0,  
      CutBasedIdMedium       = 1UL<< 1,  
      CutBasedIdMediumPrompt = 1UL<< 2,  // medium with IP cuts
      CutBasedIdTight        = 1UL<< 3,  
      CutBasedIdGlobalHighPt = 1UL<< 4,  // high pt muon for Z',W' (better momentum resolution)
      CutBasedIdTrkHighPt    = 1UL<< 5,  // high pt muon for boosted Z (better efficiency)
      PFIsoVeryLoose         = 1UL<< 6,  // reliso<0.40
      PFIsoLoose             = 1UL<< 7,  // reliso<0.25
      PFIsoMedium            = 1UL<< 8,  // reliso<0.20
      PFIsoTight             = 1UL<< 9,  // reliso<0.15
      PFIsoVeryTight         = 1UL<<10,  // reliso<0.10
      TkIsoLoose             = 1UL<<11,  // reliso<0.10
      TkIsoTight             = 1UL<<12,  // reliso<0.05
      SoftCutBasedId         = 1UL<<13,  
      SoftMvaId              = 1UL<<14,  
      MvaLoose               = 1UL<<15,  
      MvaMedium              = 1UL<<16,  
      MvaTight               = 1UL<<17,
      MiniIsoLoose           = 1UL<<18,  // reliso<0.40
      MiniIsoMedium          = 1UL<<19,  // reliso<0.20
      MiniIsoTight           = 1UL<<20,  // reliso<0.10
      MiniIsoVeryTight       = 1UL<<21   // reliso<0.05

    };
    
    bool passed( unsigned int selection ) const { return (selectors_ & selection)==selection; }
    bool passed( Selector selection ) const { return passed(static_cast<unsigned int>(selection)); }
    unsigned int selectors() const { return selectors_; }
    void setSelectors( unsigned int selectors ){ selectors_ = selectors; }
    void setSelector(Selector selector, bool passed){ 
      if (passed) selectors_ |= selector;
      else selectors_ &= ~selector;
    }
    
  private:

    bool isGlobalMuon_;
    bool isSoftMuon_;
    bool isTrackerMuon_;
    
    float normalizedChi2_;
    float ipsig_;

    int numberOfValidHits_;
    int numberOfValidMuonHits_;
    int numberOfMatchedStations_;
    int numberOfValidPixelHits_;
    int trackerLayersWithMeasurement_;

    unsigned int selectors_;
  };
}

#endif

