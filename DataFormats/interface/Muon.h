#ifndef CATTools_Muon_H
#define CATTools_Muon_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

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

    void setNormalizedChi2(float d) { normalizedChi2_ = d; }
    void setNumberOfValidHits(int i) { numberOfValidHits_ = i; }
    void setNumberOfValidMuonHits(int i) { numberOfValidMuonHits_ = i; }
    void setNumberOfMatchedStations(int i) { numberOfMatchedStations_ = i; }
    void setNumberOfValidPixelHits(int i) { numberOfValidPixelHits_ = i; }
    void setTackerLayersWithMeasurement(int i) { trackerLayersWithMeasurement_ = i; }

    void setIpSignficance(float ipsig) {ipsig_ = ipsig;}

    float scaleFactor(const std::string& name, int sign = 0) const;
    
    

  private:

    bool isGlobalMuon_;
    bool isSoftMuon_;

    float normalizedChi2_;
    float ipsig_;

    int numberOfValidHits_;
    int numberOfValidMuonHits_;
    int numberOfMatchedStations_;
    int numberOfValidPixelHits_;
    int trackerLayersWithMeasurement_;
  };
}

#endif

