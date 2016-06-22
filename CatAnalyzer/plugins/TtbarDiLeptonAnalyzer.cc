#include "CATTools/CatAnalyzer/interface/dileptonCommon.h"

using namespace std;
using namespace cat;

class TtbarDiLeptonAnalyzer : public dileptonCommon
{
public:
  explicit TtbarDiLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarDiLeptonAnalyzer() {};

private:
};

TtbarDiLeptonAnalyzer::TtbarDiLeptonAnalyzer(const edm::ParameterSet& iConfig):
  dileptonCommon(iConfig)
{
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TtbarDiLeptonAnalyzer);
