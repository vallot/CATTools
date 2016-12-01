#include "CATTools/DataFormats/interface/SecVertex.h"

using namespace cat;

/// default constructor
SecVertex::SecVertex() {
  lxy_ = -9; l3D_ = -9 ; vProb_ = -9 ; dca_ = -9 ; dca2_ = -9; dca3_= -9; cxPtHypot_ = -9; cxPtAbs_ = -9;
  leptonID1_ = -1; leptonID2_ = -1;
  trackQuality1_ = -1; trackQuality2_ =-1 ; jetDR_ = 999.; legDR_ = 999.; diffMass_ = -9;
  isMCMatch_ = false;
}

SecVertex::SecVertex(reco::VertexCompositeCandidate & aSecVertex) : reco::VertexCompositeCandidate( aSecVertex ) {
  SecVertex();
}
SecVertex::SecVertex(SecVertex const & aSecVertex) : reco::VertexCompositeCandidate( aSecVertex ) {
  lxy_ = aSecVertex.lxy(); 
  l3D_ = aSecVertex.l3D(); 
  vProb_ = aSecVertex.vProb(); 
  dca_ = aSecVertex.dca(0); 
  dca2_ = aSecVertex.dca(1); 
  dca3_ = aSecVertex.dca(2);
  cxPtHypot_ = aSecVertex.cxPtHypot();
  cxPtAbs_ = aSecVertex.cxPtAbs();
  leptonID1_ = aSecVertex.leptonID1();
  leptonID2_ = aSecVertex.leptonID2();
  trackQuality1_ = aSecVertex.trackQuality1();
  trackQuality2_ = aSecVertex.trackQuality2();
  jetDR_ = aSecVertex.JetDR();
  legDR_ = aSecVertex.LegDR();
  isMCMatch_ = aSecVertex.isMCMatch();
  diffMass_ = aSecVertex.DiffMass();
}

/// destructor
SecVertex::~SecVertex() {
}
