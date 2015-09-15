#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
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
  typedef std::vector< std::string > strings;
  typedef std::vector<std::pair<std::string, int> > stringintpairs;

  edm::EDGetTokenT<stringintpairs> catTriggerToken_;
  strings triggersToMatch_;

};

CATTriggerPacker::CATTriggerPacker(const edm::ParameterSet& pset):
  catTriggerToken_(consumes<stringintpairs>(pset.getParameter<edm::InputTag>("src"))),
  triggersToMatch_(pset.getParameter<strings>("triggersToMatch"))
{
  produces<int>("and");
  produces<int>("or");
}

void CATTriggerPacker::produce(edm::Event& event, const edm::EventSetup&)
{
  int minPS = INT_MAX, resultByAnd = 1;
  edm::Handle<stringintpairs> catTriggerHandle;
  event.getByToken(catTriggerToken_, catTriggerHandle);
  for ( const auto& catTrigger : *catTriggerHandle )
  {
    const auto& pathName = catTrigger.first;
    bool isMatched = false;
    for ( const auto& toMatch : triggersToMatch_ )
    {
      if ( pathName.find(toMatch) != std::string::npos )
      {
        isMatched = true;
        break;
      }
    }
    if ( !isMatched ) continue;
    const int result = catTrigger.second;

    resultByAnd *= result;
    if ( result ) minPS = std::min(minPS, result);
  }

  if ( minPS == INT_MAX ) minPS = 0;

  event.put(std::auto_ptr<int>(new int(resultByAnd)), "and");
  event.put(std::auto_ptr<int>(new int(minPS)), "or");

}

DEFINE_FWK_MODULE(CATTriggerPacker);
