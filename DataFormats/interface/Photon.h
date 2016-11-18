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
    
    bool isPF() const{ return isPF_; }
    bool isTight() const { return isTight_; }
    bool isMedium() const { return isMedium_; }
    bool isLoose() const { return isLoose_; }
    
    bool passMVA() const  { return passMVA_; }
    
    bool PassElectronVeto() const { return passElVeto_;}
    bool HasPixelSeed() const { return hasPixSeed_;}
    
    float SigmaiEtaiEta() const {return iEtaiEta_;}
    float r9() const {return r9_;}
    float HoverE() const { return HoverE_;}
    float SCEta() const { return scEta_;}
    float SCPhi() const { return scPhi_;}
    float SCRawEnergy() const { return scRawE_;}
    float SCPreShowerEnergy() const { return scPreShE_;}

    float smearedScale() const { return smearedScale_; }

    bool mcMatched() const { return mcMatched_; }

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
    void setIsPF(bool d) { isPF_ = d ; }
    void setIsTight(bool d) { isTight_ = d; }
    void setIsMedium(bool d) { isMedium_ = d; }
    void setIsLoose(bool d) { isLoose_ = d; }
    void setPassMVA(bool d) { passMVA_ = d; }
    void setPassElectronVeto(bool pass){ passElVeto_ = pass;}
    void setHasPixelSeed(bool hasPixSeed){hasPixSeed_ = hasPixSeed;}

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

    void setMCMatched(bool m) { mcMatched_ = m; }
    
  private:


    std::vector<pat::Photon::IdPair> photonIDs_;

    bool isPF_;
    bool isTight_;
    bool isMedium_;
    bool isLoose_;
    bool passMVA_;

    bool passElVeto_;
    bool hasPixSeed_;
    
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

    bool mcMatched_;
    
  };
}

#endif
