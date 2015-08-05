#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include <boost/regex.hpp>

#include <memory>
#include <vector>
#include <string>

class CATTriggerProducer : public edm::EDProducer
{
public:
  CATTriggerProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;
  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup) override;

private:
  typedef std::vector<double> doubles;
  typedef std::vector< std::string > strings;
  typedef std::vector< std::pair < std::string, std::string> > pairstrings;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  std::string processName_;
  pairstrings hltNames_;
  pairstrings pshltNames_;
  //HLTConfigProvider hltConfig_;
};

CATTriggerProducer::CATTriggerProducer(const edm::ParameterSet& pset)
{
  edm::InputTag hltLabel = pset.getParameter<edm::InputTag>("triggerResults");
  processName_ = hltLabel.process();
  triggerBits_ = consumes<edm::TriggerResults>(hltLabel);
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"));

  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider

  // for unprescaled triggers
  for ( auto& hltPath : pset.getParameter<strings>("unPreScaled") ){
    hltPath = boost::regex_replace(hltPath, matchVersion, "");
    std::string hltSavedAs = hltPath;
    std::replace( hltSavedAs.begin(), hltSavedAs.end(), '_', '-' );
    produces<bool>( hltSavedAs );
    hltNames_.push_back(std::make_pair(hltPath, hltSavedAs));
  }
  
  // for prescaled triggers
  for ( auto& hltPath : pset.getParameter<strings>("PreScaled") ){
    hltPath = boost::regex_replace(hltPath, matchVersion, "");
    std::string hltSavedAs = hltPath;
    std::replace( hltSavedAs.begin(), hltSavedAs.end(), '_', '-' );
    produces<int>( hltSavedAs );
    pshltNames_.push_back(std::make_pair(hltPath, hltSavedAs));
  }
}

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

void CATTriggerProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  event.getByToken(triggerBits_, triggerBits);
  event.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &trigNames = event.triggerNames(*triggerBits);

  for ( auto& hltPath : hltNames_ ){
    const strings hltPathsWithV = HLTConfigProvider::restoreVersion(trigNames.triggerNames(), hltPath.first);
    bool passTrigger = false;
    if ( hltPathsWithV.empty() ) continue;
    const std::string& trigName = hltPathsWithV[0];
    unsigned int trigIndex = trigNames.triggerIndex(trigName);
    if ( trigIndex < triggerBits->size() ){
      if ( triggerBits->accept(trigIndex) ) {
	passTrigger = true;
	break;
      }
    }
    event.put(std::auto_ptr<bool>(new bool (passTrigger)), hltPath.second);
  }

  for ( auto& hltPath : pshltNames_ ){
    const strings hltPathsWithV = HLTConfigProvider::restoreVersion(trigNames.triggerNames(), hltPath.first);
    int psValue = 0;
    if ( hltPathsWithV.empty() ) continue;
    const std::string& trigName = hltPathsWithV[0];
    unsigned int trigIndex = trigNames.triggerIndex(trigName);
    if ( trigIndex < triggerBits->size() ){
      if ( triggerBits->accept(trigIndex) ) {
	psValue = triggerPrescales->getPrescaleForIndex(trigIndex);
	break;
      }
    }
    event.put(std::auto_ptr<int>(new int (psValue)), hltPath.second);
  }
  
  // // for full list of trigger names that pass
  // for (unsigned int i = 4, n = triggerBits->size(); i < n-3; ++i) {
  //   if (triggerBits->accept(i))
  //     std::cout << i << " trigname "<<trigNames.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << std::endl;
  // }

}

DEFINE_FWK_MODULE(CATTriggerProducer);

