#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <boost/regex.hpp>
#include <vector>
#include <string>

class CATTriggerPacker : public edm::stream::EDProducer<>
{
public:
  CATTriggerPacker(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup&) override;

private:
  //typedef std::vector<bool> vbool;
  //typedef std::vector<int> vint;
  typedef std::vector<std::string> strings;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> prescaleToken_;
  strings triggersToMatch_;
  const bool combineByOr_;

};

CATTriggerPacker::CATTriggerPacker(const edm::ParameterSet& pset):
  triggerToken_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerResults"))),
  prescaleToken_(consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("triggerPrescales"))),
  triggersToMatch_(pset.getParameter<strings>("triggersToMatch")),
  combineByOr_(pset.getParameter<bool>("combineByOr"))
{
  produces<int>();
}

void CATTriggerPacker::produce(edm::Event& event, const edm::EventSetup&)
{
  using namespace std;

  int result = combineByOr_ ? -1 : 1;

  edm::Handle<edm::TriggerResults> triggerHandle;
  event.getByToken(triggerToken_, triggerHandle);
  edm::Handle<pat::PackedTriggerPrescales> prescaleHandle;
  event.getByToken(prescaleToken_, prescaleHandle);

  const auto& trigNames = event.triggerNames(*triggerHandle);

  // Keep trigger indices and ps factors
  int minPS = INT_MAX;
  std::vector<std::pair<int, int> > triggerInfo;
  for ( const auto& trigName : trigNames.triggerNames() )
  {
    const unsigned int index = trigNames.triggerIndex(trigName);
    for ( const auto& trigPattern : triggersToMatch_ )
    {
      if ( trigName.find(trigPattern) != string::npos )
      {
        const int ps = prescaleHandle->getPrescaleForIndex(index);
        triggerInfo.push_back(make_pair(index, ps));
        minPS = std::min(ps, minPS);
        break;
      }
    }
  }
  if ( minPS == INT_MAX ) minPS = 0; // This may happen if trigger names are wrong

  // Get the trigger results and combine them
  for ( const auto ti : triggerInfo )
  {
    const int index = ti.first, ps = ti.second;
    const bool isFired = triggerHandle->accept(index);

    if ( isFired )
    {
      if ( combineByOr_ ) { result = minPS; break; }
      else { result *= ps; }
    }
    else
    {
      if ( !combineByOr_ ) { result = 0; break; }
    }
  }

  if ( result == -1 ) result = 0;
  // Result should be:
  //   combine by or, anything fired => min(ps1, ps2, ...)
  //   combine by or, nothing fired  => -1 => replace to 0
  //   combine by and, anything not fired => 0
  //   combine by and, everything fired => ps1*ps2*...

  event.put(std::auto_ptr<int>(new int(result)));
}

DEFINE_FWK_MODULE(CATTriggerPacker);
