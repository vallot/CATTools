#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
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

#include "CATTools/DataFormats/interface/TriggerBits.h"

#include <boost/regex.hpp>

#include <memory>
#include <vector>
#include <string>

using namespace cat;
using namespace std;

class CATTriggerProducer : public edm::stream::EDProducer<>
{
public:
  CATTriggerProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;
  //void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup) override;

private:
  typedef std::vector<double> doubles;
  typedef std::vector< std::string > strings;
  typedef std::vector< std::pair < std::string, std::string> > pairstrings;

  std::vector<edm::EDGetTokenT<edm::TriggerResults>> flagResTokens_;
  std::map<std::string, edm::EDGetTokenT<bool>> flagBoolTokens_;
  strings flagNamesToSelect_;

  std::vector<edm::EDGetTokenT<edm::TriggerResults>> trigResTokens_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> trigPSToken_;
  strings trigNamesToSelect_;
  pairstrings hltNames_;

  //HLTConfigProvider hltConfig_;
};

CATTriggerProducer::CATTriggerProducer(const edm::ParameterSet& pset)
{
  auto hltPSet = pset.getParameter<edm::ParameterSet>("HLT");
  for ( auto x : hltPSet.getParameter<std::vector<edm::InputTag>>("triggerResults") ) {
    trigResTokens_.push_back(consumes<edm::TriggerResults>(x));
  }
  trigObjToken_ = consumes<pat::TriggerObjectStandAloneCollection>(hltPSet.getParameter<edm::InputTag>("objects"));
  trigPSToken_ = consumes<pat::PackedTriggerPrescales>(hltPSet.getParameter<edm::InputTag>("prescales"));
  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider
  for ( auto& name : hltPSet.getParameter<strings>("select") ) {
    trigNamesToSelect_.push_back(name);
    string hltPath = boost::regex_replace(name, matchVersion, "");
    std::string hltSavedAs = hltPath;
    hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());
    std::cout << " " << hltPath << std::endl;
    hltNames_.push_back(std::make_pair(hltPath, hltSavedAs));
  }

  auto flagPSet = pset.getParameter<edm::ParameterSet>("Flag");
  flagNamesToSelect_ = flagPSet.getParameter<strings>("names");
  for ( auto x : flagPSet.getParameter<std::vector<edm::InputTag>>("Flags") ) {
    flagResTokens_.push_back(consumes<edm::TriggerResults>(x));
  }
  for ( auto x : flagPSet.getParameter<std::vector<edm::InputTag>>("bools") ) {
    flagBoolTokens_.insert(std::make_pair("Flag_"+x.label(), consumes<bool>(x)));
    flagNamesToSelect_.push_back("Flag_"+x.label());
  }

  produces<cat::TriggerBits>();
  produces<pat::TriggerObjectStandAloneCollection >();
}

/*
void CATTriggerProducer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  // bool changed = true;
  // if ( !hltConfig_.init(run, eventSetup, processName_, changed) )
  // {
  //   edm::LogError("CATTriggerProducer") << "HLT config extraction failure with process name " << processName_;
  // }

  // if ( changed )
  // {
  //   edm::LogError("CATTriggerProducer") << "HLT config has changed " << processName_;
  //  // The HLT config has actually changed wrt the previous Run, hence rebook your
  //  // histograms or do anything else dependent on the revised HLT config
  // }
}
*/

void CATTriggerProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::unique_ptr<cat::TriggerBits> out_trigBits(new cat::TriggerBits());
  std::unique_ptr<pat::TriggerObjectStandAloneCollection> out_objects(new pat::TriggerObjectStandAloneCollection());

  std::map<std::string, int> results;

  // Open RECO filter flag handles
  std::map<string, edm::Handle<bool>> flagBoolHandles;
  for ( auto key=flagBoolTokens_.begin(); key!=flagBoolTokens_.end(); ++key ) {
    edm::Handle<bool> handle;
    if ( !event.getByToken(key->second, handle) ) continue;
    flagBoolHandles.insert(std::make_pair(key->first, handle));
  };
  std::vector<edm::Handle<edm::TriggerResults>> flagResHandles;
  for ( auto token : flagResTokens_ ) {
    edm::Handle<edm::TriggerResults> handle;
    if ( !event.getByToken(token, handle) ) continue;
    flagResHandles.push_back(handle);
  }

  for ( string flagName : flagNamesToSelect_ ) {
    int result = 0;
    
    auto flagBool = flagBoolHandles.find(flagName);
    if ( flagBool != flagBoolHandles.end() ) {
      result = *(flagBool->second);
    }

    for ( auto handle : flagResHandles ) {
      const auto& filterNames = event.triggerNames(*handle);
      unsigned int index = filterNames.triggerIndex(flagName);
      if ( index < filterNames.size() and handle->accept(index) ) {
        result = 1;
        break;
      }
    }

    results[flagName] = result;
  }

  // Trigger results and objects
  for ( auto trigResToken : trigResTokens_ ) {
    edm::Handle<edm::TriggerResults> trigResHandle;
    if ( !event.getByToken(trigResToken, trigResHandle) ) continue;
    const edm::TriggerNames &trigNames = event.triggerNames(*trigResHandle);

    // saving trigger info as preScale value (int)
    edm::Handle<pat::PackedTriggerPrescales> trigPSHandle;
    event.getByToken(trigPSToken_, trigPSHandle);
    for ( auto& hltPath : hltNames_ ){
      const strings hltPathsWithV = HLTConfigProvider::restoreVersion(trigNames.triggerNames(), hltPath.first);
      int psValue = 0;
      if ( hltPathsWithV.empty() ){
        //std::cout << "Warning:: trigger does not exist "<< hltPath.first << std::endl;
        continue;
      }
      const std::string& trigName = hltPathsWithV[0];
      unsigned int trigIndex = trigNames.triggerIndex(trigName);
      if ( trigIndex < trigResHandle->size() ){
        if ( trigResHandle->accept(trigIndex) ) {
          psValue = trigPSHandle->getPrescaleForIndex(trigIndex);
        }
      }
      results[hltPath.second] = psValue;
    }

    edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjHandle;
    event.getByToken(trigObjToken_, trigObjHandle);

    for ( pat::TriggerObjectStandAlone trigObj : *trigObjHandle ) { // note: not "const &" since we want to call unpackPathNames
      bool keepTriggerObject = false;
      trigObj.unpackPathNames(trigNames);
      std::vector<std::string> pathNamesAll  = trigObj.pathNames(false);
      for ( unsigned int h = 0, n = pathNamesAll.size(); h < n; ++h ) {
        for ( auto& trigName : trigNamesToSelect_ ){
          if ( pathNamesAll[h].find(trigName) == 0 and
               trigObj.hasPathName( pathNamesAll[h], true, true ) ) {
              keepTriggerObject = true;
              break;
          }
        }
      }
      if ( keepTriggerObject ) {
        trigObj.packPathNames(trigNames);
        out_objects->push_back(trigObj);
      }
    }
  }

  out_trigBits->set(results);
  event.put(std::move(out_trigBits));
  event.put(std::move(out_objects));
}

DEFINE_FWK_MODULE(CATTriggerProducer);

