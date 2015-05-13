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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include <boost/regex.hpp>

#include <memory>
#include <vector>
#include <string>

class RecoEventInfoProducer : public edm::EDProducer
{
public:
  RecoEventInfoProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;
  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup) override;

private:
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  std::string processName_;
  std::map<std::string, strings> hltGroup_;
  //HLTConfigProvider hltConfig_;

};

RecoEventInfoProducer::RecoEventInfoProducer(const edm::ParameterSet& pset)
{
  vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));
  edm::InputTag hltLabel = pset.getParameter<edm::InputTag>("triggerResults");
  processName_ = hltLabel.process();
  triggerBits_ = consumes<edm::TriggerResults>(hltLabel);
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"));

  produces<int>("pvN");
  produces<double>("pvX");
  produces<double>("pvY");
  produces<double>("pvZ");

  edm::ParameterSet hltSet = pset.getParameter<edm::ParameterSet>("HLT");
  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider
  for ( auto& hltSetName : hltSet.getParameterNamesForType<strings>() )
    {
      const std::string hltGroupName = hltSetName;
      strings& hltPaths = hltGroup_[hltGroupName];
      hltPaths = hltSet.getParameter<strings>(hltGroupName);
      for ( auto& hltPath : hltPaths )
	{
	  hltPath = boost::regex_replace(hltPath, matchVersion, "");
	}

      produces<int>("HLT"+hltSetName);
    }
}

void RecoEventInfoProducer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  // bool changed = true;
  // if ( !hltConfig_.init(run, eventSetup, processName_, changed) ) 
  // {
  //   edm::LogError("RecoEventInfoProducer") << "HLT config extraction failure with process name " << processName_;
  // }

  // if ( changed ) 
  // {
  //   edm::LogError("RecoEventInfoProducer") << "HLT config has changed " << processName_;
  //  // The HLT config has actually changed wrt the previous Run, hence rebook your
  //  // histograms or do anything else dependent on the revised HLT config
  // }
}

void RecoEventInfoProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);

  const int nPV = vertexHandle->size();
  double pvX = 0, pvY = 0, pvZ = 0;
  if ( nPV > 0 )
    {
      pvX = vertexHandle->at(0).x();
      pvY = vertexHandle->at(0).y();
      pvZ = vertexHandle->at(0).z();
    }
  event.put(std::auto_ptr<int>(new int(nPV)), "pvN");
  event.put(std::auto_ptr<double>(new double(pvX)), "pvX");
  event.put(std::auto_ptr<double>(new double(pvY)), "pvY");
  event.put(std::auto_ptr<double>(new double(pvZ)), "pvZ");

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  event.getByToken(triggerBits_, triggerBits);
  event.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &trigNames = event.triggerNames(*triggerBits);

  for ( auto key = hltGroup_.begin(); key != hltGroup_.end(); ++key )
    {
      const std::string& hltGroupName = key->first;
      const strings& hltPaths = key->second;

      int psValue = 0;
      for ( auto& hltPath : hltPaths )
	{
	  const strings hltPathsWithV = HLTConfigProvider::restoreVersion(trigNames.triggerNames(), hltPath);
	  if ( hltPathsWithV.empty() ) continue;
	  const std::string& trigName = hltPathsWithV[0];
	  unsigned int trigIndex = trigNames.triggerIndex(trigName);
	  if ( trigIndex < triggerBits->size() ){
	    if ( triggerBits->accept(trigIndex) ) {
	      // const std::pair<int,int> prescales = hltConfig_.prescaleValues(event, eventSetup, trigName);
	      // psValue = prescales.first * prescales.second;
	      psValue = triggerPrescales->getPrescaleForIndex(trigIndex);
	      break;
	    }
	  }
	}
      event.put(std::auto_ptr<int>(new int (psValue)), "HLT"+hltGroupName);
    }

  // // for full list of trigger names that pass
  // for (unsigned int i = 4, n = triggerBits->size(); i < n-3; ++i) {
  //   if (triggerBits->accept(i))
  //     std::cout << i << " trigname "<<trigNames.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << std::endl;
  // }

}

DEFINE_FWK_MODULE(RecoEventInfoProducer);

