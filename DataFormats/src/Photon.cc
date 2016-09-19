#include "CATTools/DataFormats/interface/Photon.h"

using namespace cat;

/// default constructor
Photon::Photon() {
}

Photon::Photon(const reco::LeafCandidate & aPhoton) :
  Particle( aPhoton ),
  isPF_(false),
  isTight_(false),
  isMedium_(false),
  isLoose_(false),
  passMVA_(false),
  passElVeto_(false),
  hasPixSeed_(false),
  rhoIso_(0),
  chargedHadronIso_(0),
  puChargedHadronIso_(0),
  neutralHadronIso_(0),
  photonIso_(0),
  iEtaiEta_(0),
  r9_(0),
  HoverE_(0),
  scEta_(0), scPhi_(0), scRawE_(0), scPreShE_(0),
  mcMatched_(false)
{}

/// destructor
Photon::~Photon() {
}



float Photon::photonID(const std::string& name) const {
  for (std::vector<pat::Photon::IdPair>::const_iterator it = photonIDs_.begin(), ed = photonIDs_.end(); it != ed; ++it) {
    if (it->first.find(name) != std::string::npos)
      return it->second;
  }
  cms::Exception ex("Key not found");
  ex << "cat::Photon: the ID " << name << " can't be found in this cat::Photon.\n";
  ex << "The available IDs are: ";
  for (std::vector<pat::Photon::IdPair>::const_iterator it = photonIDs_.begin(), ed = photonIDs_.end(); it != ed; ++it) {
    ex << "'" << it->first << "' ";
  }
  ex << ".\n";
  throw ex;
}

float Photon::chargedHadronIsoWithEA(){
  return std::max( (float)0.0, chargedHadronIso_ - rhoIso_*getEffArea(charged, abs(scEta_)));
}

float Photon::neutralHadronIsoWithEA(){
  return std::max( (float)0.0, neutralHadronIso_ - rhoIso_*getEffArea(neutral, abs(scEta_)));
}

float Photon::photonIsoWithEA(){
  return std::max( (float)0.0, photonIso_ - rhoIso_*getEffArea(photon, abs(scEta_)));
}

float Photon::getEffArea(IsoType iso_type, float scEta) 
{

  /// Taken from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2                                                                                                                                                                            

  //iso_type = 0 (charged hadron) ; 1 (neutral hadron) ; 2 (photons)                                                                                                                                                                                                  
  float absEta = std::abs(scEta);

  if(iso_type == charged) return 0.;
  else if (iso_type == neutral){
    if ( 0.0000 >= absEta && absEta < 1.0000 ) return 0.0599;
    if ( 1.0000 >= absEta && absEta < 1.4790 ) return 0.0819;
    if ( 1.4790 >= absEta && absEta < 2.0000 ) return 0.0696;
    if ( 2.0000 >= absEta && absEta < 2.2000 ) return 0.0360;
    if ( 2.2000 >= absEta && absEta < 2.3000 ) return 0.0360;
    if ( 2.3000 >= absEta && absEta < 2.4000 ) return 0.0462;
    if ( 2.4000 >= absEta && absEta < 5.0000 ) return 0.0656;
  }
  else if (iso_type == photon){
    if ( 0.0000 >= absEta && absEta < 1.0000 ) return 0.1271;
    if ( 1.0000 >= absEta && absEta < 1.4790 ) return 0.1101;
    if ( 1.4790 >= absEta && absEta < 2.0000 ) return 0.0756;
    if ( 2.0000 >= absEta && absEta < 2.2000 ) return 0.1175;
    if ( 2.2000 >= absEta && absEta < 2.3000 ) return 0.1498;
    if ( 2.3000 >= absEta && absEta < 2.4000 ) return 0.1857;
    if ( 2.4000 >= absEta && absEta < 5.0000 ) return 0.2183;
  }

  return 0;
}


