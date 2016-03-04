#ifndef CATTools_GenWeights_H
#define CATTools_GenWeights_H

#include <vector>
#include <string>

// Define typedefs for convenience
namespace cat {
  class GenWeights;
}

namespace cat {

class GenWeightInfo {
public:
  GenWeightInfo() {};
  virtual ~GenWeightInfo() {};

  int nGroups() const { return typeNames_.size(); }
  std::string typeName(int i) const { return typeNames_.at(i); }
  std::string combineMethod(int i) const { return combineMethods_.at(i); }
  std::vector<std::string> weightParams(int i) const { return weightParams_.at(i); }

  int print() const;

  void addWeightGroup(const std::string typeName, const std::string combineBy, const std::vector<std::string> weightParams, const std::vector<std::string> weightKey);

  enum KnownTypes { NONE=-1, Nominal, PDF, ScaleUp, ScaleDown };
  static KnownTypes toKnownType(std::string typeName);

private:
  std::vector<std::string> typeNames_;
  std::vector<std::string> combineMethods_;
  std::vector<std::vector<std::string> > weightKeys_;
  std::vector<std::vector<std::string> > weightParams_;

};

class GenWeights {
public:
  GenWeights();
  virtual ~GenWeights() {};

  // Getters
  float genWeight() const { return genWeight_; }
  float lheWeight() const { return lheWeight_; }
  std::vector<float> scaleUpWeights() const { return scaleUpWeights_; }
  std::vector<float> scaleDownWeights() const { return scaleDownWeights_; }
  std::vector<float> pdfWeights() const { return pdfWeights_; }
  std::vector<float> otherWeights() const { return otherWeights_; }

  int id1() const { return id1_; }
  int id2() const { return id2_; }
  float x1() const { return x1_; }
  float x2() const { return x2_; }
  float qScale() const { return qScale_; }

  // Setters
  void setInfo(const int id1, const int id2, const float x1, const float x2, const float qScale);

  void setGenWeight(const float w) { genWeight_ = w; }
  void setLHEWeight(const float w) { lheWeight_ = w; }
  void addScaleUpWeight(const float w) { scaleUpWeights_.push_back(w); }
  void addScaleDownWeight(const float w) { scaleDownWeights_.push_back(w); }
  void addPDFWeight(const float w) { pdfWeights_.push_back(w); }
  void addOtherWeight(const float w) { otherWeights_.push_back(w); }

private:
  int id1_, id2_;
  float x1_, x2_, qScale_;

  float genWeight_, lheWeight_;
  std::vector<float> scaleUpWeights_, scaleDownWeights_;
  std::vector<float> pdfWeights_;
  std::vector<float> otherWeights_;

};

}

#endif
