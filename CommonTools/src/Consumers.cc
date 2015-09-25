#include "CATTools/CommonTools/interface/Consumers.h"

namespace cat {
  // Build dummy instances here
  FlatConsumers<bool> boolCSet;
  FlatConsumers<int> intCSet;
  FlatConsumers<double> doubleCSet;
  FlatConsumers<float> floatCSet;
  FlatConsumers<std::string> stringCSet;
  
  VectorConsumers<bool> vboolCSet;
  VectorConsumers<int> vintCSet;
  VectorConsumers<double> vdoubleCSet;
  VectorConsumers<float> vfloatCSet;
  VectorConsumers<std::string> vstringCSet;
}
