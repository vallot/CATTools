#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
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
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_, secondaryTriggerToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> prescaleToken_;
  strings triggersToMatch_;
  bool combineByOr_;

};

CATTriggerBitCombiner::CATTriggerBitCombiner(const edm::ParameterSet& pset):
  doFilter_(pset.getParameter<bool>("doFilter")),
  triggerToken_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerResults"))),
  prescaleToken_(consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("triggerPrescales"))),
  triggersToMatch_(pset.getParameter<strings>("triggersToMatch"))
{
  if ( pset.existsAs<edm::InputTag>("secondaryTriggerResults") )
  {
    secondaryTriggerToken_ = consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("secondaryTriggerResults"));
  }

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

  int result = combineByOr_ ? -1 : 1;

  edm::Handle<edm::TriggerResults> triggerHandle;
  if ( !event.getByToken(triggerToken_, triggerHandle) )
  {
    event.getByToken(secondaryTriggerToken_, triggerHandle);
  }
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

  if ( !doFilter_ ) return true;
  return (result != 0);
}

DEFINE_FWK_MODULE(CATTriggerBitCombiner);
