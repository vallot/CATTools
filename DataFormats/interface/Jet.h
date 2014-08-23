#ifndef CATTools_Jet_H
#define CATTools_Jet_H 

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <string>
#include <boost/array.hpp>

// Define typedefs for convenience
namespace cat {
  class Jet;
  typedef std::vector<Jet>              JetCollection;
  typedef edm::Ref<JetCollection>       JetRef;
  typedef edm::RefVector<JetCollection> JetRefVector;
}

namespace cat {

  class Jet : public Particle{
  public:
    Jet();
    Jet(const reco::LeafCandidate & aJet); 
    virtual ~Jet();

    bool LooseId() const { return LooseId_; }

    /// \return btag discriminator
    double btag(unsigned int index = 0) const{ return index < btag_.size() ? btag_.at(index) : -9.9; }
    double btag(const char* s) const;
    double bDiscriminator(const char* s) const{ return btag(s); }

    enum BTagWP { TCHPT, JPL, JPM, JPT, CSVL, CSVM, CSVT,  CSVV1L, CSVV1M, CSVV1T,  CSVSLV1L, CSVSLV1M, CSVSLV1T, CSVIVFV2L, CSVIVFV2M, CSVIVFV2T  };
    bool btagWP(BTagWP wp) const ;
    bool btagWP(const std::string &wp) const ;
    bool btagWP(const char *wp) const ;

    /// \return secondary vertex b-tagging information
    Float_t secvtxMass() const { return secvtxMass_ ; }
    Float_t Lxy() const { return Lxy_ ; }
    Float_t LxyErr() const { return LxyErr_; }

    /// \return the matched MC parton flavour (from the shower, used e.g. for b-tagging)
    Int_t partonFlavour() const{ return partonFlavour_;}
    /// \return the pdgId of the matched MC parton from hard scattering (i.e. the closest quark or gluon of status == 3)
    Int_t partonPdgId() const{ return partonPdgId_;}

    void setbTag( const int & i, const double & d, const std::string & name ) { 
      btag_[i] = d;
      btagNames_[i] = name;
    }

    void setSecVtxMass(Float_t f) { secvtxMass_ = f;}
    void setLxy(Float_t f) { Lxy_ = f;}
    void setLxyErr(Float_t f) { LxyErr_ = f;}
    void setPartonFlavour(Int_t i) { partonFlavour_ = i; }
    void setPartonPdgId(Int_t i) { partonPdgId_ = i; }
    void setLooseId(bool id) { LooseId_ = id; } 

  private:

    bool LooseId_; 
 
    /// b tagging discriminators
    typedef boost::array<double,16> TagArray;
    typedef boost::array<std::string,TagArray::static_size> TagNameArray;
    TagArray btag_;
    TagNameArray btagNames_;

    /// b tagging information
    Float_t secvtxMass_;
    Float_t Lxy_;
    Float_t LxyErr_;

    //parton flavour
    Int_t partonFlavour_;
    Int_t partonPdgId_;

  };
}

#endif
