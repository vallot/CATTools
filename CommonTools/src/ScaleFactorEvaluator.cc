#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include <cassert>

using namespace cat;

void ScaleFactorEvaluator::set(const std::vector<double>& xbins,
                               const std::vector<double>& ybins,
                               const std::vector<double>& values,
                               const std::vector<double>& errors)
{
  xbins_ = xbins;
  ybins_ = ybins;
  values_ = values;
  errors_ = errors;

  const unsigned int n = (xbins.size()-1)*(ybins.size()-1);
  // FIXME : check that these bins are monolothically increasing
  assert(values.size() == n);
  assert(errors.size() == n);

  // For cache
  width_ = xbins_.size()-1;
}

double ScaleFactorEvaluator::operator()(const double x, const double y, const double shift) const
{
  // Filter out UF and OF
  if ( x < (*xbins_.begin()) or x >= (*(xbins_.end()-1)) ) return 1;
  if ( y < (*ybins_.begin()) or y >= (*(ybins_.end()-1)) ) return 1;

  // Special note: std::lower_bound DOES NOT GIVE lower bound of the bin,
  // it gives the first element which satiesfiles (x < VECTOR_ELEMENT)
  // - as it is desigend to properly work on integers
  auto xbin = std::lower_bound(xbins_.begin(), xbins_.end(), x);
  auto ybin = std::lower_bound(ybins_.begin(), ybins_.end(), y);
  if ( (*xbin) != x ) --xbin;
  if ( (*ybin) != y ) --ybin;

  const int column = xbin-xbins_.begin();
  const int row = ybin-ybins_.begin();

  const int bin = row*width_+column;
  const double value = values_.at(bin);
  const double error = errors_.at(bin);

  return std::max(0.0, value+shift*error);
}

double ScaleFactorEvaluator::getScaleFactor(const cat::Particle& p, const int pid, const double shift) const
{
  const int aid = abs(p.pdgId());
  if ( aid == pid ) {
    const double x = p.pt(), y = p.eta();
    
    auto xbin = std::lower_bound(xbins_.begin(), xbins_.end(), x);
    if ( xbin == xbins_.end() or xbin+1 == xbins_.end() ) return 1;
    auto ybin = std::lower_bound(ybins_.begin(), ybins_.end(), y);
    if ( ybin == ybins_.end() or ybin+1 == ybins_.end() ) return 1;

    const int column = xbin-xbins_.begin();
    const int row = ybin-ybins_.begin();

    const int bin = row*width_+column;
    const double value = values_.at(bin);
    const double error = errors_.at(bin);

    return std::max(0.0, value+shift*error);
  }
  return 1;
}
