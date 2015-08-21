#include "CATTools/DataFormats/interface/MET.h"

using namespace cat;

MET::MET(const reco::LeafCandidate & aMet) : Particle( aMet ) {
}

MET::MET(const reco::LeafCandidate & aMet, float sumEt) : Particle( aMet ), sumEt_(sumEt){
}
