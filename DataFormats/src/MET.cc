#include "CATTools/DataFormats/interface/MET.h"

using namespace cat;

/// default constructor
MET::MET() {
}

MET::MET(float px, float py, float sumEt) {
  px_ = px; py_ = py; sumEt_ = sumEt;
}

/// destructor
MET::~MET() {
}
