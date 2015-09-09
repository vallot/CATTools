#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

class CATTriggerPacker : public edm::stream::EDProducer<>
{
public:
  CATTriggerPacker(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup&) override;

private:
  typedef std::vector<bool> vbool;
  typedef std::vector<int> vint;

  std::vector<edm::EDGetTokenT<bool> > boolTokens_;
  std::vector<edm::EDGetTokenT<int> > intTokens_;

};

CATTriggerPacker::CATTriggerPacker(const edm::ParameterSet& pset)
{
  for ( auto x : pset.getParameter<std::vector<edm::InputTag> >("srcs") )
  {
    boolTokens_.push_back(consumes<bool>(x));
    intTokens_.push_back(consumes<int>(x));
  }

  produces<int>("and");
  produces<int>("or");
}

void CATTriggerPacker::produce(edm::Event& event, const edm::EventSetup&)
{
  int minPS = INT_MAX, resultByAnd = 1;
  for ( const auto& tok : boolTokens_ )
  {
    edm::Handle<bool> handle;
    if ( !event.getByToken(tok, handle) ) continue;
    resultByAnd *= *handle;
    if ( *handle ) minPS = 1;
  }

  for ( const auto& tok : intTokens_ )
  {
    edm::Handle<int> handle;
    if ( !event.getByToken(tok, handle) ) continue;
    resultByAnd *= *handle;
    if ( *handle ) minPS = std::min(minPS, *handle);
  }

  if ( minPS == INT_MAX ) minPS = 0;

  event.put(std::auto_ptr<int>(new int(resultByAnd)), "and");
  event.put(std::auto_ptr<int>(new int(minPS)), "or");
}

DEFINE_FWK_MODULE(CATTriggerPacker);
