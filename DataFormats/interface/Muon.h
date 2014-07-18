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

    double chargedHadronIso() const { return chargedHadronIso_;  }
    double puChargedHadronIso() const { return puChargedHadronIso_; }
    double neutralHadronIso() const { return neutralHadronIso_; }
    double photonIso() const { return photonIso_; }

    bool isTightMuon() const { return isTightMuon_; }
    bool isLooseMuon() const { return isLooseMuon_; } 
    bool isSoftMuon() const { return isSoftMuon_; } 
 
    void setChargedHadronIso(double i) { chargedHadronIso_ = i; }
    void setPUChargedHadronIso(double i) { puChargedHadronIso_ = i; }
    void setNeutralHadronIso(double i) { neutralHadronIso_ = i; }
    void setPhotonIso(double i) { photonIso_ = i; }

    void setIsTightMuon(bool d) { isTightMuon_ = d; }
    void setIsLooseMuon(bool d) { isLooseMuon_ = d; }
    void setIsSoftMuon(bool d) { isSoftMuon_ = d; }
     
 
  private:

    double chargedHadronIso_;
    double puChargedHadronIso_;
    double neutralHadronIso_;
    double photonIso_;

    bool isTightMuon_; 
    bool isLooseMuon_; 
    bool isSoftMuon_;

  };
}

#endif
