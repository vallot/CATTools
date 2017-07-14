#ifndef CATTools_CommonTools_GenParticleHelper_H
#define CATTools_CommonTools_GenParticleHelper_H

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TLorentzVector.h"
#include <vector>

/*
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/LorentzVector.h"
*/

namespace cat {

std::vector<reco::GenParticleRef> selectGenParticles(const edm::Handle<reco::GenParticleCollection>& srcColl,
                                                     const unsigned int pdgId, const std::vector<unsigned int>& motherIds);

bool isMatchedByDRDPt(const reco::Candidate::LorentzVector& a, const reco::Candidate::LorentzVector& b, bool exact);

bool isLast(const reco::GenParticle& p);

}

#endif
