#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
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

class LumiMaskProducer : public edm::stream::EDProducer<>
{
public:
  LumiMaskProducer(const edm::ParameterSet& pset);
  ~LumiMaskProducer() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  std::vector<edm::LuminosityBlockRange> luminositySectionsBlockRanges_;
};

LumiMaskProducer::LumiMaskProducer(const edm::ParameterSet& pset)
{
  luminositySectionsBlockRanges_ = pset.getUntrackedParameter< std::vector<edm::LuminosityBlockRange> >("LumiSections");
  produces<int>("");
}

void LumiMaskProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<int> passLumiSelection(new int(0));
  if ( event.isRealData() ){
    edm::RunNumber_t             kRun  = event.id().run();
    edm::LuminosityBlockNumber_t kLumi = event.id().luminosityBlock();
    for(std::vector<edm::LuminosityBlockRange>::iterator oneLumiRange = luminositySectionsBlockRanges_.begin();
  	oneLumiRange != luminositySectionsBlockRanges_.end(); ++oneLumiRange){
      if ( kRun < (*oneLumiRange).startRun() ) {
  	break;
      }
      if ( (*oneLumiRange).endRun() < kRun ) continue;
      if ( (*oneLumiRange).startLumi() <= kLumi && kLumi <= (*oneLumiRange).endLumi() ){
  	*passLumiSelection = 1;
  	break;
      }
    }
  }
  event.put(passLumiSelection, "");
}

DEFINE_FWK_MODULE(LumiMaskProducer);
