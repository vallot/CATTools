#include "CATTools/DataFormats/interface/MET.h"

using namespace cat;

MET::MET(const reco::LeafCandidate & aMet) :
  Particle( aMet ),
  sumEt_(0),
  rawMET_(0),
  unclusteredEnUp_px_(0),unclusteredEnUp_py_(0),unclusteredEnUp_sumet_(0),
  unclusteredEnDown_px_(0),unclusteredEnDown_py_(0),unclusteredEnDown_sumet_(0),
  jetEnUp_px_(0),jetEnUp_py_(0),jetEnUp_sumet_(0),
  jetEnDown_px_(0),jetEnDown_py_(0),jetEnDown_sumet_(0),
  jetResUp_px_(0),jetResUp_py_(0),jetResUp_sumet_(0),
  jetResDown_px_(0),jetResDown_py_(0),jetResDown_sumet_(0)
{}

MET::MET(const reco::LeafCandidate & aMet, float sumEt) :
  Particle( aMet ),
  sumEt_(sumEt),
  rawMET_(0),
  unclusteredEnUp_px_(0),unclusteredEnUp_py_(0),unclusteredEnUp_sumet_(0),
  unclusteredEnDown_px_(0),unclusteredEnDown_py_(0),unclusteredEnDown_sumet_(0),
  jetEnUp_px_(0),jetEnUp_py_(0),jetEnUp_sumet_(0),
  jetEnDown_px_(0),jetEnDown_py_(0),jetEnDown_sumet_(0),
  jetResUp_px_(0),jetResUp_py_(0),jetResUp_sumet_(0),
  jetResDown_px_(0),jetResDown_py_(0),jetResDown_sumet_(0)
{}
