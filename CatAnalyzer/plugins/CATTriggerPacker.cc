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
  typedef std::vector<bool> vbool;
  typedef std::vector<int> vint;
  typedef std::vector< std::string > strings;
  typedef std::pair<std::string, int> stringintpair;
  typedef std::vector< stringintpair > stringintpairs;

  edm::EDGetTokenT<stringintpairs > InputTriggerLabel_;
  std::vector<std::string>  *triggerNames_;

  std::vector<edm::EDGetTokenT<bool> > boolTokens_;
  std::vector<edm::EDGetTokenT<int> > intTokens_;

};

CATTriggerPacker::CATTriggerPacker(const edm::ParameterSet& pset):
InputTriggerLabel_(consumes< stringintpairs >(pset.getParameter<edm::InputTag>("InputTriggerLabel")))
{
  for ( auto x : pset.getParameter<std::vector<edm::InputTag> >("srcs") )
  {
    boolTokens_.push_back(consumes<bool>(x));
    intTokens_.push_back(consumes<int>(x));
  }

  produces<int>("and");
  produces<int>("or");

  triggerNames_           = new std::vector<std::string>;
  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider

  //std::cout << "List of Triggers to Save" << std::endl;
  for ( auto& hltPath : pset.getParameter<strings>("hltPathNames") ){
    hltPath = boost::regex_replace(hltPath, matchVersion, "");
    std::string hltSavedAs = hltPath;
    hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());
    //std::cout << " " << hltPath << std::endl;
    produces<int >( hltSavedAs );
    triggerNames_->push_back(hltSavedAs);
  }

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


  edm::Handle<stringintpairs> Triggers;
  event.getByToken(InputTriggerLabel_, Triggers);
  unsigned int n = Triggers->size();

  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider
  for ( unsigned int ip=0; ip<n; ++ip ) 
  {
    const stringintpair& p = (*Triggers)[ip];
    std::string name_   = p.first;
    int         isPath_ = p.second;
    //std::cout << name_ <<":"<< isPath_ << std::endl;
    name_ = boost::regex_replace(name_, matchVersion, "");
    std::string hltSavedAs = name_;
    hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());

    for(unsigned int it2=0;it2<triggerNames_->size();it2++)
    {
      std::size_t found2 = triggerNames_->at(it2).find(hltSavedAs);
      if (found2!=std::string::npos && hltSavedAs.length()==triggerNames_->at(it2).length() )
          event.put(std::auto_ptr<int>(new int (isPath_)), hltSavedAs);    
          //std::cout << hltSavedAs <<":"<< isPath_ << std::endl;
    }
  }

}

DEFINE_FWK_MODULE(CATTriggerPacker);
