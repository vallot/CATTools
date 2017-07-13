#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "CATTools/DataFormats/interface/Trigger.h"

#include <boost/regex.hpp>

#include <memory>
#include <vector>
#include <string>

class CATTriggerProducer : public edm::one::EDProducer<edm::BeginLuminosityBlockProducer>
{
public:
  CATTriggerProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup&) override;
  void beginLuminosityBlockProduce(edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup) override;

private:
  typedef std::vector<double> doubles;
  typedef std::vector< std::string > strings;
  typedef std::vector< std::pair < std::string, std::string> > pairstrings;
  typedef pat::TriggerObjectStandAloneCollection TrigObjColl;

  std::vector<edm::EDGetTokenT<edm::TriggerResults>> trigResTokens_;
  edm::EDGetTokenT<TrigObjColl> trigObjToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> trigPSToken_;
  std::vector<std::string> trigPrefixes_;
  std::vector<edm::EDGetTokenT<edm::TriggerResults>> flagTokens_;
  std::vector<std::string> flagNames_;
  std::map<std::string, edm::EDGetTokenT<bool>> flagBools_;
  std::vector<std::string> trigProcs_;

  HLTConfigProvider hltConfig_;
  cat::TriggerResValues filterBits_;
};

CATTriggerProducer::CATTriggerProducer(const edm::ParameterSet& pset)
{
  auto flagPSet = pset.getParameter<edm::ParameterSet>("flags");
  flagNames_ = flagPSet.getParameter<strings>("names");
  for ( auto& x : flagPSet.getParameter<std::vector<edm::InputTag>>("triggerResults") ) {
    flagTokens_.push_back(consumes<edm::TriggerResults>(x));
  }
  for ( auto& x : flagPSet.getParameter<std::vector<edm::InputTag>>("bools") ) {
    if ( x.instance().empty() ) flagBools_[x.label()] = consumes<bool>(x);
    else flagBools_[x.label()+"_"+x.instance()] = consumes<bool>(x);
  }

  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider

  auto trigPSet = pset.getParameter<edm::ParameterSet>("trigger");
  trigObjToken_ = consumes<TrigObjColl>(trigPSet.getParameter<edm::InputTag>("objects"));
  trigPSToken_ = consumes<pat::PackedTriggerPrescales>(trigPSet.getParameter<edm::InputTag>("prescales"));
  for ( auto x : trigPSet.getParameter<std::vector<edm::InputTag>>("triggerResults") ) {
    trigResTokens_.push_back(consumes<edm::TriggerResults>(x));
    trigProcs_.push_back(x.process());
  }
  trigPrefixes_ = trigPSet.getParameter<strings>("prefix");

/*
  for ( auto& hltPath : pset.getParameter<strings>("hltPathNames") ){
    hltPath = boost::regex_replace(hltPath, matchVersion, "");
    std::string hltSavedAs = hltPath;
    hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());
    std::cout << " " << hltPath << std::endl;
    produces<int >( hltSavedAs );
    hltNames_.push_back(std::make_pair(hltPath, hltSavedAs));
  }
*/

  produces<cat::TriggerNames, edm::InLumi>();
  produces<cat::TriggerBits>();
  produces<TrigObjColl>();
}

using namespace std;
void CATTriggerProducer::beginLuminosityBlockProduce(edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
  filterBits_.clear();

  for ( std::string name : flagNames_ ) {
    if ( name.find("Flag_") != 0 ) name = "Flag_"+name;
    filterBits_[name] = 0;
  }
  for ( auto key : flagBools_ ) {
    std::string name = key.first;
    if ( name.find("Flag_") != 0 ) name = "Flag_"+name;
    filterBits_[name] = 0;
  }

  const auto& run = lumi.getRun();
  for ( auto procName : trigProcs_ ) {
    bool isChanged = true;
    if ( hltConfig_.init(run, eventSetup, procName, isChanged) ) break;
  }
  for ( auto name : hltConfig_.triggerNames() ) {
    bool skipPath = true;
    for ( auto& prefix : trigPrefixes_ ) {
      if ( name.find(prefix) == 0 ) {
        skipPath = false;
        break;
      }
    }
    if ( skipPath ) continue;

    filterBits_[name] = 0;
  }

  lumi.put(std::auto_ptr<cat::TriggerNames>(new cat::TriggerNames(filterBits_)));
}

void CATTriggerProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  // Initialize all filter results to 0
  for ( auto key = filterBits_.begin(); key != filterBits_.end(); ++key ) {
    key->second = 0;
  }

  // Collect event filter results from TriggerResults
  // Give the priority to the last in
  for ( auto token : flagTokens_ ) {
    edm::Handle<edm::TriggerResults> flagHandle;
    if ( !event.getByToken(token, flagHandle) ) continue;

    auto filterNames = event.triggerNames(*flagHandle);
    for ( auto flagName : flagNames_ ) {
      unsigned int trigIndex = filterNames.triggerIndex(flagName);
      if ( trigIndex >= flagHandle->size() ) continue;

      if ( flagName.find("Flag_") != 0 ) flagName = "Flag_"+flagName;
      if ( filterBits_.find(flagName) == filterBits_.end() ) continue;
      filterBits_[flagName] = flagHandle->accept(trigIndex);
    }
  }

  // Collect event filter results from EDProducers
  // Overrides bits from TriggerResults
  for ( auto key : flagBools_ ) {
    std::string name = key.first;
    edm::Handle<bool> flagHandle;
    if ( !event.getByToken(key.second, flagHandle) ) continue;

    if ( name.find("Flag_") != 0 ) name = "Flag_"+name;
    if ( filterBits_.find(name) == filterBits_.end() ) continue;
    filterBits_[name] = *flagHandle;
  }

  // Load HLT flags. Use the first successful one.
  edm::Handle<edm::TriggerResults> trigResHandle;
  for ( auto token : trigResTokens_ ) {
    if ( event.getByToken(token, trigResHandle) ) break;
  }
  auto& trigNames = event.triggerNames(*trigResHandle);
  edm::Handle<pat::PackedTriggerPrescales> trigPSHandle;
  event.getByToken(trigPSToken_, trigPSHandle);

  for ( auto pathName : trigNames.triggerNames() ) {
    unsigned int trigIndex = trigNames.triggerIndex(pathName);
    if ( trigIndex >= trigResHandle->size() ) continue;

    bool skipPath = true;
    for ( auto& prefix : trigPrefixes_ ) {
      if ( pathName.find(prefix) == 0 ) {
        skipPath = false;
        break;
      }
    }
    if ( skipPath ) continue;

    int psValue = 0;
    if ( trigResHandle->accept(trigIndex) ) {
      psValue = trigPSHandle->getPrescaleForIndex(trigIndex);
    }

    if ( filterBits_.find(pathName) == filterBits_.end() ) continue;

    filterBits_[pathName] = std::min(USHRT_MAX, psValue);
  }

  // Load trigger objects
  edm::Handle<TrigObjColl> trigObjHandle;
  event.getByToken(trigObjToken_, trigObjHandle);
  std::unique_ptr<TrigObjColl> out_trigObjs(new TrigObjColl());
  for ( pat::TriggerObjectStandAlone trigObj : *trigObjHandle ) { // note: not "const &" since we want to call unpackPathNames
    trigObj.unpackPathNames(trigNames);
    bool keepObject = false;
    for ( auto pathName : trigObj.pathNames(false) ) {
      for ( auto& prefix : trigPrefixes_ ) {
        if ( pathName.find(prefix) == 0 and
             trigObj.hasPathName(pathName, true, true) ) {
          keepObject = true;
          break;
        }
      }
      if ( keepObject ) break;
    }
    if ( keepObject ) {
      trigObj.packPathNames(trigNames);
      out_trigObjs->push_back(trigObj);
    }
  }

  event.put(std::move(std::unique_ptr<cat::TriggerBits>(new cat::TriggerBits(filterBits_))));
  event.put(std::move(out_trigObjs));
}

DEFINE_FWK_MODULE(CATTriggerProducer);

