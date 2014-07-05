#include "../interface/GenTop.h"

using namespace cat;

bool GenTop::isHadronicWellSeparated(float deltaMin_) const
{
	if(!isHadronic_) return false;
	float d1 = ROOT::Math::VectorUtil::DeltaR(quark_,quarkBar_);
	float d2 = ROOT::Math::VectorUtil::DeltaR(quark_,bquark_);
	float d3 = ROOT::Math::VectorUtil::DeltaR(bquark_,quarkBar_);
	if(d1 > deltaMin_ && d2 > deltaMin_ && d3 > deltaMin_) return true;
	return false;
}
  
float GenTop::DeltaRMinHadronicTop() const
{
	if(!isHadronic_) return -999.;
	float dmin = -999;
	float d1 = ROOT::Math::VectorUtil::DeltaR(quark_,quarkBar_);
	float d2 = ROOT::Math::VectorUtil::DeltaR(quark_,bquark_);
	float d3 = ROOT::Math::VectorUtil::DeltaR(bquark_,quarkBar_);
	if(d1<d2) dmin = d1;
	else dmin = d2;
	if(d3<dmin) dmin = d3;
	return dmin;
}
 
ClassImp(GenTop)
 
