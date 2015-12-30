#include <vector>
#include <utility>
#include <algorithm>

namespace cat {

class ScaleFactorGetter
{
public:
  void set(const std::vector<double>& xbins,
           const std::vector<double>& ybins,
           const std::vector<double>& values,
           const std::vector<double>& errors);
  double operator()(const double x, const double y, const double shift = 0) const;
  
private:
  std::vector<double> xbins_, ybins_;
  std::vector<double> values_;
  std::vector<double> errors_;

  int width_;
};

}
