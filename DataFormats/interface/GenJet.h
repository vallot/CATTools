#ifndef CATTools_GenJet_H
#define CATTools_GenJet_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "CATTools/DataFormats/interface/MCParticle.h"

#include <string>
#include <boost/array.hpp>

// Define typedefs for convenience
namespace cat {
  class GenJet;
  typedef std::vector<GenJet>              GenJetCollection;
  typedef edm::Ref<GenJetCollection>       GenJetRef;
  typedef edm::RefVector<GenJetCollection> GenJetRefVector;
}

namespace cat {

  class GenJet : public Particle{
  public:
    GenJet();
    GenJet(const reco::GenJet & aGenJet);
    virtual ~GenJet();

    const MCParticle hadron() const {return hadron_; }
    void setHadron( MCParticle Had ) { hadron_ = Had; }	
    /// \return the matched MC parton flavour (from the shower, used e.g. for b-tagging)
    int partonFlavour() const{ return partonFlavour_;}
    /// \return the pdgId of the matched MC parton from hard scattering (i.e. the closest quark or gluon of status == 3)
    int partonPdgId() const{ return partonPdgId_;}
    void setPartonFlavour(int i) { partonFlavour_ = i; }
    void setPartonPdgId(int i) { partonPdgId_ = i; }

  private:

    MCParticle hadron_;
    //parton flavour
    int partonFlavour_;
    int partonPdgId_;

  };
}

#endif
