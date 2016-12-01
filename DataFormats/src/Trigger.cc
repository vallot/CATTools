#include "CATTools/DataFormats/interface/Trigger.h"

#include <iostream>
#include <algorithm>

using namespace cat;
using namespace std;

TriggerNames::TriggerNames()
{
}

void TriggerNames::set(const cat::TriggerResValues& results)
{
  names_.clear();
  for ( auto key = results.begin(); key != results.end(); ++key ) {
    names_.push_back(key->first);
  }
}

int TriggerNames::index(const std::string& name) const
{
  auto match = std::equal_range(names_.begin(), names_.end(), name).first;
  if ( match == names_.end() ) return -1;

  return match-names_.begin();
}

int TriggerNames::print() const
{
  for ( auto& name : names_ ) cout << name << endl;
  return 0;
}

TriggerBits::TriggerBits()
{
}

void TriggerBits::set(const cat::TriggerResValues& results)
{
  for ( auto key = results.begin(); key != results.end(); ++key ) {
    values_.push_back(key->second);
  }
}

