#ifndef HLTAnalyzer_h
#define HLTAnalyzer_h

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "CATTools/DataFormats/interface/Run.h"
#include "CATTools/DataFormats/interface/Event.h"

#include <iomanip>

using namespace cat;

class HLTAnalyzer
{
	
public:
	
	HLTAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig) :
		verbosity_(0)
		,triggerResultsTag1st_(producersNames.getParameter<edm::InputTag> ("hltProducer1st"))
		,triggerResultsTag2nd_(producersNames.getParameter<edm::InputTag> ("hltProducer2nd"))
		,triggerResultsTag3rd_(producersNames.getParameter<edm::InputTag> ("hltProducer3rd"))
		,triggerResultsTag4th_(producersNames.getParameter<edm::InputTag> ("hltProducer4th"))
  		,triggerNames_()
		,doHLT_(myConfig.getUntrackedParameter<bool>("doHLT",false))
  		,nEvents_(0)
  		,nWasRun_(0)
  		,nAccept_(0)
  		,nErrors_(0)
  		,hltWasRun_(0)
  		,hltAccept_(0)
  		,hltErrors_(0)
  		,hltNames_(0)
		{;}
	
	~HLTAnalyzer() {;}
	
	void setVerbosity(int verbosity) {verbosity_ = verbosity; };
	void init(const edm::Event& iEvent, Event* rootEvent);
	void process(const edm::Event& iEvent, Event* rootEvent);
	void printStats();
	void copySummary(Run* runInfos);
	
/*	unsigned int  nHLTPaths() const { return hltNames_.size(); }

	unsigned int  nEvents() const { return nEvents_; }
	unsigned int  nWasRun() const { return nWasRun_; }
	unsigned int  nAccept() const { return nAccept_; }
	unsigned int  nErrors() const { return nErrors_; }

	std::vector<unsigned int> hltWasRun() const { return hltWasRun_; }
	std::vector<unsigned int> hltAccept() const { return hltAccept_; }
	std::vector<unsigned int> hltErrors() const { return hltErrors_; }
	std::vector<std::string>  hltNames() const { return hltNames_; }

	unsigned int hltWasRun(unsigned ipath) const { return (hltWasRun_.size()>ipath ?  hltWasRun_.at(ipath) : 0 ); }
	unsigned int hltAccept(unsigned ipath) const { return (hltAccept_.size()>ipath ?  hltAccept_.at(ipath) : 0 ); }
	unsigned int hltErrors(unsigned ipath) const { return (hltErrors_.size()>ipath ?  hltErrors_.at(ipath) : 0 ); }
	std::string hltNames(unsigned ipath) const { return (hltNames_.size()>ipath ?  hltNames_.at(ipath) : "noname" ); }
*/

private:
	int verbosity_;
	
	edm::InputTag triggerResultsTag_;		// Input tag for TriggerResults, final choice
	edm::InputTag triggerResultsTag1st_;	// Input tag for TriggerResults, 1st choice
	edm::InputTag triggerResultsTag2nd_;	// Input tag for TriggerResults, 2nd choice
	edm::InputTag triggerResultsTag3rd_;	// Input tag for TriggerResults, 3rd choice
	edm::InputTag triggerResultsTag4th_;	// Input tag for TriggerResults, 4th choice

	edm::TriggerNames triggerNames_;			// TriggerNames class

	bool doHLT_;

	unsigned int  nEvents_;								// number of events processed

	unsigned int  nWasRun_;								// # where at least one HLT was run
	unsigned int  nAccept_;								// # of accepted events
	unsigned int  nErrors_;								// # where at least one HLT had error

	std::vector<unsigned int> hltWasRun_;			// # where HLT[i] was run
	std::vector<unsigned int> hltAccept_;			// # of events accepted by HLT[i]
	std::vector<unsigned int> hltErrors_;			// # of events with error in HLT[i]
	std::vector<std::string>  hltNames_;			// name of each HLT algorithm


	// new HLTInfo container

	vector<cat::HLTInfo> hltInfos_;
	
};

#endif
