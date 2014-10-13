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

    const MCParticle BHadron() const {return BHadron_; }
    const MCParticle CHadron() const {return CHadron_; }

    void setBHadron( MCParticle BHad ) { BHadron_ = BHad; }	
    void setCHadron( MCParticle CHad ) { CHadron_ = CHad; }	

  private:

    MCParticle BHadron_;
    MCParticle CHadron_;

  };
}

#endif
