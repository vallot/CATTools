#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include <vector>

namespace reco {
  typedef std::vector<CompositePtrCandidate> CompositePtrCandidateCollection;
  typedef edm::Ref<CompositePtrCandidateCollection> CompositePtrCandidateRef;
  typedef edm::RefVector<CompositePtrCandidateCollection> CompositePtrCandidateRefVector;
}

namespace {

struct CATTools_PhysicsAnalysis_DataFormats {
  edm::Wrapper<reco::CompositePtrCandidateCollection> a1;
  edm::reftobase::Holder<reco::Candidate, reco::CompositePtrCandidateRef> a2;
  edm::reftobase::RefHolder<reco::CompositePtrCandidateRef> a3;
  edm::reftobase::VectorHolder<reco::Candidate, reco::CompositePtrCandidateRefVector> a4;
  edm::reftobase::RefVectorHolder<reco::CompositePtrCandidateRefVector> a5;

};

}
