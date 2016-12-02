#include "FWCore/Framework/interface/stream/EDFilter.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CATTools/DataFormats/interface/Trigger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <vector>
#include <string>

class CATTriggerBitCombiner : public edm::stream::EDFilter<>
{
public:
  CATTriggerBitCombiner(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;

private:
  //typedef std::vector<bool> vbool;
  //typedef std::vector<int> vint;
  typedef std::vector<std::string> strings;
  const bool doFilter_;
  edm::EDGetTokenT<cat::TriggerNames> trigNamesToken_;
  edm::EDGetTokenT<cat::TriggerBits> trigBitsToken_;
  strings triggersToMatch_;
  bool combineByOr_;

};

CATTriggerBitCombiner::CATTriggerBitCombiner(const edm::ParameterSet& pset):
  doFilter_(pset.getParameter<bool>("doFilter")),
  trigNamesToken_(consumes<cat::TriggerNames, edm::InLumi>(pset.getParameter<edm::InputTag>("src"))),
  trigBitsToken_(consumes<cat::TriggerBits>(pset.getParameter<edm::InputTag>("src"))),
  triggersToMatch_(pset.getParameter<strings>("triggersToMatch"))
{
  auto combineBy = pset.getParameter<std::string>("combineBy");
  std::transform(combineBy.begin(), combineBy.end(), combineBy.begin(), ::toupper);
  if ( combineBy == "OR" ) combineByOr_ = true;
  else if ( combineBy == "AND" ) combineByOr_ = false;
  else edm::LogError("CATTriggerBitCombiner") << "Wrong input to \"combinedBy\" option, it was \"" << pset.getParameter<std::string>("combineBy") << ".\n"
                                              << "This should be chosen among (\"and\", \"or\")\n";
  produces<int>();
}

bool CATTriggerBitCombiner::filter(edm::Event& event, const edm::EventSetup&)
{
  using namespace std;

  edm::Handle<cat::TriggerNames> trigNamesHandle;
  event.getLuminosityBlock().getByToken(trigNamesToken_, trigNamesHandle);

  edm::Handle<cat::TriggerBits> trigBitsHandle;
  event.getByToken(trigBitsToken_, trigBitsHandle);

  // Keep trigger indices and ps factors
  std::vector<int> results;
  bool hasFailed = false;
  const auto trigNames = trigNamesHandle->names();
  for ( int index=0, nTrig=trigNames.size(); index<nTrig; ++index ) {
    const auto& trigName = trigNames[index];
    const int accPS = trigBitsHandle->result(index);
    for ( const auto& trigPattern : triggersToMatch_ ) {
      if ( trigName.find(trigPattern) != 0 ) continue;
      if ( accPS == 0 ) hasFailed = true;
      else results.push_back(accPS);
    }
  }

  int result = 0;
  if ( combineByOr_ ) {
    if ( results.empty() ) result = 0;
    else result = *(std::min_element(results.begin(), results.end()));
  }
  else {
    if ( results.empty() or hasFailed ) result = 0;
    else {
      result = std::accumulate(results.begin(), results.end(), 1, std::multiplies<int>());
    }
  }

  // Result should be:
  //   combine by or, anything fired => min(ps1, ps2, ...)
  //   combine by or, nothing fired  => 0
  //   combine by and, anything not fired => 0
  //   combine by and, everything fired => ps1*ps2*...

  event.put(std::auto_ptr<int>(new int(result)));

  if ( !doFilter_ ) return true;
  return (result != 0);
}

DEFINE_FWK_MODULE(CATTriggerBitCombiner);
