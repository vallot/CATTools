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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include <boost/regex.hpp>

#include <memory>
#include <vector>
#include <string>

class CATTriggerProducer : public edm::stream::EDProducer<>
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
  edm::EDGetTokenT<edm::TriggerResults> metFilterBitsPAT_;
  edm::EDGetTokenT<edm::TriggerResults> metFilterBitsRECO_;

  pairstrings hltNames_;
  pairstrings metFilterNames_;
  //HLTConfigProvider hltConfig_;
};

CATTriggerProducer::CATTriggerProducer(const edm::ParameterSet& pset):
  triggerBits_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("bits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"))),
  metFilterBitsPAT_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("metFilterBitsPAT"))),
  metFilterBitsRECO_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("metFilterBitsRECO")))
{
  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider

  std::cout << "List of Triggers to Save" << std::endl;
  for ( auto& hltPath : pset.getParameter<strings>("hltPathNames") ){
    hltPath = boost::regex_replace(hltPath, matchVersion, "");
    std::string hltSavedAs = hltPath;
    hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());
    std::cout << " " << hltPath << std::endl;
    produces<int >( hltSavedAs );
    hltNames_.push_back(std::make_pair(hltPath, hltSavedAs));
  }
  produces<std::vector< std::pair < std::string, int > >>();

  for ( auto& hltPath : pset.getParameter<strings>("metFilterNames") ){
    std::cout << " " << hltPath << std::endl;
    produces<bool >( hltPath );
    metFilterNames_.push_back(std::make_pair("Flag_"+hltPath, hltPath));
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
    int psValue = 0;
    if ( hltPathsWithV.empty() ){
      //std::cout << "Warning:: trigger does not exist "<< hltPath.first << std::endl;
      continue;
    }
    const std::string& trigName = hltPathsWithV[0];
    unsigned int trigIndex = trigNames.triggerIndex(trigName);
    if ( trigIndex < triggerBits->size() ){
      if ( triggerBits->accept(trigIndex) ) {
	psValue = triggerPrescales->getPrescaleForIndex(trigIndex);
      }
    }
    event.put(std::auto_ptr<int>(new int (psValue)), hltPath.second);    
  }

  // save only ele and mu triggers that pass
  std::vector< std::pair < std::string, int > > *alltriggers = new std::vector< std::pair < std::string, int > >();
  for( unsigned int i=0; i<trigNames.size(); ++i ){
    if (trigNames.triggerName(i).find("HLT_Ele") == 0 
	|| trigNames.triggerName(i).find("HLT_DoubleEle") == 0 
	|| trigNames.triggerName(i).find("HLT_IsoMu") == 0 
	|| trigNames.triggerName(i).find("HLT_Mu") == 0 ){
      if ( triggerBits->accept(i) ) {
	int psValue = int(triggerBits->accept(i)) * triggerPrescales->getPrescaleForIndex(i);
	alltriggers->push_back(std::make_pair(trigNames.triggerName(i), psValue));
      }
    }
  }
  event.put(std::auto_ptr<std::vector< std::pair < std::string, int > >>(alltriggers));

  // save filter info
  edm::Handle<edm::TriggerResults> metFilterBits;
  if (!event.getByToken(metFilterBitsPAT_, metFilterBits)){
    event.getByToken(metFilterBitsRECO_, metFilterBits);
  }
  
  const edm::TriggerNames &metFilterNames = event.triggerNames(*metFilterBits);

  for ( auto& hltPath : metFilterNames_ ){
    bool passMet = false;
    unsigned int trigIndex = metFilterNames.triggerIndex(hltPath.first);
    if ( trigIndex < metFilterBits->size() ){
      if ( metFilterBits->accept(trigIndex) )
	passMet = true;
    }
    event.put(std::auto_ptr<bool>(new bool (passMet)), hltPath.second);
  }
  
  // // for full list of metFilterger names that pass
  // for (unsigned int i = 0, n = metFilterBits->size(); i < n-3; ++i) {
  //   std::cout << i << " metFiltername "<<metFilterNames.triggerName(i)<< std::endl;
  // }

}

DEFINE_FWK_MODULE(CATTriggerProducer);

