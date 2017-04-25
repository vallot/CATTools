#ifndef CATTools_Photon_H
#define CATTools_Photon_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

// Define typedefs for convenience
namespace cat {
  class Photon;
  typedef std::vector<Photon>              PhotonCollection;
  typedef edm::Ref<PhotonCollection>       PhotonRef;
  typedef edm::RefVector<PhotonCollection> PhotonRefVector;
}

namespace cat {

  class Photon : public Particle{
  public:
    Photon();
    Photon(const reco::LeafCandidate & aPhoton);
    virtual ~Photon();


    enum IsoType         {charged,
                          neutral,
                          photon};


    float photonID(const std::string& name) const;
    float photonID(const char* name) const { return photonID( std::string(name) );} 


    float chargedHadronIso() const { return chargedHadronIso_;  }
    float puChargedHadronIso() const { return puChargedHadronIso_; }
    float neutralHadronIso() const { return neutralHadronIso_; }
    float photonIso() const { return photonIso_; }
    float rhoIso()  const { return rhoIso_;}
    float chargedHadronIsoWithEA(); 
    float neutralHadronIsoWithEA(); 
    float photonIsoWithEA(); 
    
    bool isPF() const{ return idBits_[0]; }
    bool isTight() const { return idBits_[1]; }
    bool isMedium() const { return idBits_[2]; }
    bool isLoose() const { return idBits_[3]; }
    
    bool passMVA() const  { return idBits_[4]; }
    
    bool PassElectronVeto() const { return idBits_[5];}
    bool HasPixelSeed() const { return idBits_[6];}
    
    float SigmaiEtaiEta() const {return iEtaiEta_;}
    float r9() const {return r9_;}
    float HoverE() const { return HoverE_;}
    float SCEta() const { return scEta_;}
    float SCPhi() const { return scPhi_;}
    float SCRawEnergy() const { return scRawE_;}
    float SCPreShowerEnergy() const { return scPreShE_;}

    float smearedScale() const { return smearedScale_; }

    bool mcMatched() const { return idBits_[7]; }

    float getEffArea(IsoType iso_type, float scEta) ;
    /// Set class variables

    void setPhotonIDs(const std::vector<pat::Photon::IdPair> & ids) { photonIDs_ = ids; }
    void setPhotonID(pat::Photon::IdPair ids) { photonIDs_.push_back(ids); }

    // Isolation variables
    void setRho(double rho) {rhoIso_ = rho;}
    void setchargedHadronIso(float ch){ chargedHadronIso_ = ch;}
    void setneutralHadronIso(float nh){ neutralHadronIso_ = nh;}
    void setphotonIso(float ph){ photonIso_ = ph;}
    void setpuhargedHadronIso(float pu) {puChargedHadronIso_ = pu;} 
    
    // ID variables
    void setIsPF(bool d) { idBits_[0] = d ; }
    void setIsTight(bool d) { idBits_[1] = d; }
    void setIsMedium(bool d) { idBits_[2] = d; }
    void setIsLoose(bool d) { idBits_[3] = d; }
    void setPassMVA(bool d) { idBits_[4] = d; }
    void setPassElectronVeto(bool pass){ idBits_[5] = pass;}
    void setHasPixelSeed(bool hasPixSeed){ idBits_[6] = hasPixSeed;}

    // Shower shape
    void setSigmaiEtaiEta(float ieteieta){ iEtaiEta_ = ieteieta;}
    void setR9(float varr9){ r9_ = varr9;}
    void setHoverE(float hovere){HoverE_ = hovere;}

    // SC variables
    void setSCEta(float eta){ scEta_ = eta;}
    void setSCPhi(float phi){ scPhi_ = phi;}
    void setSCrawEnergy(float raw_e){ scRawE_ = raw_e;}
    void setSCPreShowerEnergy(float pre_E){ scPreShE_ = pre_E;}

    void setSmearedScale(const float scale) { smearedScale_ = scale; }

    void setMCMatched(bool m) { idBits_[7] = m; }
    
  private:


    std::vector<pat::Photon::IdPair> photonIDs_;

    std::bitset<8> idBits_; // isPF, isTight, isMedium, isLoose, passMVA, passElVeto, hasPixSeed, mcMatched
    
    float rhoIso_;
    float chargedHadronIso_;
    float puChargedHadronIso_;
    float neutralHadronIso_;
    float photonIso_;
    
    float iEtaiEta_;
    float r9_;
    float HoverE_;
    float scEta_, scPhi_, scRawE_, scPreShE_;

    float smearedScale_;
  };
}

#endif
