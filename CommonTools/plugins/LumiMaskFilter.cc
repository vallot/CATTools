#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class LumiMaskFilter : public edm::stream::EDFilter<>
{
public:
  LumiMaskFilter(const edm::ParameterSet& pset);
  ~LumiMaskFilter() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  std::vector<edm::LuminosityBlockRange> luminositySectionsBlockRanges_;
  const bool doFilter_, acceptOnFail_;
};

LumiMaskFilter::LumiMaskFilter(const edm::ParameterSet& pset):
  doFilter_(pset.getUntrackedParameter<bool>("doFilter", false)), // Default behaviour is to use it as a producer
  acceptOnFail_(pset.getUntrackedParameter<bool>("acceptOnFail", false)) // Default behaviour is to accept event if it is in the range
{
  luminositySectionsBlockRanges_ = pset.getUntrackedParameter< std::vector<edm::LuminosityBlockRange> >("LumiSections");
  produces<int>("");
}

bool LumiMaskFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<int> passLumiSelection(new int(0));
  if ( event.isRealData() ) {
    edm::RunNumber_t             kRun  = event.id().run();
    edm::LuminosityBlockNumber_t kLumi = event.id().luminosityBlock();
    for ( auto oneLumiRange : luminositySectionsBlockRanges_ ) {
      if ( kRun < oneLumiRange.startRun() ) break;
      if ( oneLumiRange.endRun() < kRun ) continue;
      if ( oneLumiRange.startLumi() <= kLumi && kLumi <= oneLumiRange.endLumi() ) {
        *passLumiSelection = 1;
        break;
      }
    }
  }
  const bool isPassing = (*passLumiSelection == 1);
  event.put(passLumiSelection, "");
  if ( doFilter_ ) {
    // XOR selection, but let us keep the logic to be simple to everybody
    if      ( isPassing and !acceptOnFail_ ) return true;
    else if ( !isPassing and acceptOnFail_ ) return true;
    else return false;
  }
  return true;
}

DEFINE_FWK_MODULE(LumiMaskFilter);
