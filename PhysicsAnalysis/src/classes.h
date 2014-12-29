#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include <vector>

namespace {

struct CATTools_PhysicsAnalysis_DataFormats {
    reco::CompositePtrCandidate a_;
    std::vector<reco::CompositePtrCandidate> av;
    edm::Wrapper<reco::CompositePtrCandidate> aw;
    edm::Wrapper<std::vector<reco::CompositePtrCandidate> > avw;
    edm::Ptr<reco::CompositePtrCandidate> aPtr;
};

}
