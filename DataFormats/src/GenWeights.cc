#include "CATTools/DataFormats/interface/GenWeights.h"
#include <algorithm>
#include <iostream>

using namespace std;
using namespace cat;

void GenWeightInfo::addWeightGroup(const string typeName, const string combineBy, const vector<string> weightParams, const vector<string> weightKeys)
{
  typeNames_.push_back(typeName);
  combineMethods_.push_back(combineBy);
  weightParams_.push_back(weightParams);
  weightKeys_.push_back(weightKeys);
}

GenWeightInfo::KnownTypes GenWeightInfo::toKnownType(string typeName)
{
  std::transform(typeName.begin(), typeName.end(), typeName.begin(), ::tolower);

  if ( typeName == "nominal" ) return Nominal;
  else if ( typeName == "pdf" ) return PDF;
  else if ( typeName == "scaleup" ) return ScaleUp;
  else if ( typeName == "scaledown" ) return ScaleDown;

  return NONE;
}

/// default constructor
GenWeights::GenWeights()
{
  lheWeight_ = genWeight_ = 1;
  id1_ = id2_ = 0;
  x1_ = x2_ = qScale_ = 0;
}

void GenWeights::setInfo(const int id1, const int id2, const float x1, const float x2, const float qScale)
{
  id1_ = id1;
  id2_ = id2;
  x1_ = x1;
  x2_ = x2;
  qScale_ = qScale;
}

int GenWeightInfo::print() const
{
  cout << "------------------------------------------------\n";
  for ( int i=0, n=typeNames_.size(); i<n; ++i ) {
    cout << typeNames_[i] << ' ' << combineMethods_[i] << endl;
    for ( int j=0, m=weightParams_[i].size(); j<m; ++j ) {
      cout << "  -> " << weightParams_[i][j] << " access using " << weightKeys_[i][j] << endl;
    }
  }
  cout << "------------------------------------------------\n";

  return 0;
}

