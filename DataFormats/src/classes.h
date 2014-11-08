#include "DataFormats/Common/interface/Wrapper.h"

#include "CATTools/DataFormats/interface/Particle.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Photon.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "CATTools/DataFormats/interface/SecVertex.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>

#include <TMatrixD.h>

namespace {
  struct CATTools_DataFormats {

    cat::Particle a_;
    std::vector<cat::Particle> av;
    edm::Wrapper<cat::Particle> aw;
    edm::Wrapper<std::vector<cat::Particle> > avw;
    edm::Ptr<cat::Particle> aPtr;

    cat::Muon m_;
    std::vector<cat::Muon> mv;
    edm::Wrapper<cat::Muon> mw;
    edm::Wrapper<std::vector<cat::Muon> > mvw;
    edm::Ptr<cat::Muon> mPtr;

    cat::Electron e_;
    std::vector<cat::Electron> ev;
    edm::Wrapper<cat::Electron> ew;
    edm::Wrapper<std::vector<cat::Electron> > evw;
    edm::Ptr<cat::Electron> ePtr;

    cat::Photon p_;
    std::vector<cat::Photon> pv;
    edm::Wrapper<cat::Photon> pw;
    edm::Wrapper<std::vector<cat::Photon> > pvw;
    edm::Ptr<cat::Photon> pPtr;

    cat::Jet j_;
    std::vector<cat::Jet> jv;
    edm::Wrapper<cat::Jet> jw;
    edm::Wrapper<std::vector<cat::Jet> > jvw;
    edm::Ptr<cat::Jet> jPtr;

    cat::MET met_;
    std::vector<cat::MET> metv;
    edm::Wrapper<cat::MET> metw;
    edm::Wrapper<std::vector<cat::MET> > metvw;
    edm::Ptr<cat::MET> metPtr;

    cat::MCParticle ma_;
    std::vector<cat::MCParticle> mav;
    edm::Wrapper<cat::MCParticle> maw;
    edm::Wrapper<std::vector<cat::MCParticle> > mavw;
    edm::Ptr<cat::MCParticle> maPtr;

    cat::GenJet gj_;
    std::vector<cat::GenJet> gjv;
    edm::Wrapper<cat::GenJet> gjw;
    edm::Wrapper<std::vector<cat::GenJet> > gjvw;
    edm::Ptr<cat::GenJet> gjPtr;

    cat::GenTop tgj_;
    std::vector<cat::GenTop> tgjv;
    edm::Wrapper<cat::GenTop> tgjw;
    edm::Wrapper<std::vector<cat::GenTop> > tgjvw;
    edm::Ptr<cat::GenTop> tgjPtr;

    cat::SecVertex sv_;
    std::vector<cat::SecVertex> svv;
    edm::Wrapper<cat::SecVertex> svw;
    edm::Wrapper<std::vector<cat::SecVertex> > svvw;
    edm::Ptr<cat::SecVertex> svPtr;

  };
  

}
