#include "../interface/HLTAnalyzer.h"

using namespace std;
using namespace cat;

int getIndexForRun ( const edm::Event& iEvent, vector<cat::HLTInfo> infos ) {
  
  for (unsigned int i=0; i<infos.size(); i++) {

    cat::HLTInfo hltInfo_ = infos[i];

    if ( iEvent.id().run() == hltInfo_.RunID() )
    
      return i;

  }

  return -1;

}

/*void HLTAnalyzer::init(const edm::Event& iEvent, cat::Event* rootEvent)
{
       
   edm::Handle<edm::TriggerResults> trigResults;
   try {iEvent.getByLabel(triggerResultsTag1st_,trigResults);} catch (...) {;}
   if (!trigResults.isValid())
     doHLT_=false;
   else
     doHLT_=true;     
}
*/

void HLTAnalyzer::process(const edm::Event& iEvent, cat::Event* rootEvent)
{
	nEvents_++;
	
	if(doHLT_)
	{
	  
	  int index = getIndexForRun(iEvent,hltInfos_);
		
//		cout << "index " << index << " and hltInfos_ " << hltInfos_.size() << endl;

	  
	  if ( index == -1 || index < (int)hltInfos_.size()-1 ) {// no info for this run yet -> create a new entry

	    edm::Handle<edm::TriggerResults> trigResults, trigResults1st, trigResults2nd, trigResults3rd, trigResults4th;
	    try {iEvent.getByLabel(triggerResultsTag1st_,trigResults1st);} catch (...) {;}
	    try {iEvent.getByLabel(triggerResultsTag2nd_,trigResults2nd);} catch (...) {;}
	    try {iEvent.getByLabel(triggerResultsTag3rd_,trigResults3rd);} catch (...) {;}
	    try {iEvent.getByLabel(triggerResultsTag4th_,trigResults4th);} catch (...) {;}
	    if(trigResults1st.isValid())
	      {
				trigResults = trigResults1st;
				triggerResultsTag_ = triggerResultsTag1st_;
	      }
	    else if(trigResults2nd.isValid())
	      {
				trigResults = trigResults2nd;
				triggerResultsTag_ = triggerResultsTag2nd_;
	      }
	    else if(trigResults3rd.isValid())
	      {
				trigResults = trigResults3rd;
				triggerResultsTag_ = triggerResultsTag3rd_;
	      }
	    else if(trigResults4th.isValid())
	    {
	      trigResults = trigResults4th;
	      triggerResultsTag_ = triggerResultsTag4th_;
			}
			
	    triggerNames_=iEvent.triggerNames(*trigResults);
	    hltNames_=triggerNames_.triggerNames();
	    const unsigned int n(hltNames_.size());
	    hltWasRun_.resize(n);
	    hltAccept_.resize(n);
	    hltErrors_.resize(n);
	    for (unsigned int i=0; i!=n; ++i)
	      {
					hltWasRun_[i]=0;
					hltAccept_[i]=0;
					hltErrors_[i]=0;
	      }
	    
	    if (index == -1) {
	      hltInfos_.push_back(cat::HLTInfo(iEvent.id().run(),hltNames_,hltWasRun_,hltAccept_,hltErrors_));
	      index = getIndexForRun(iEvent,hltInfos_); // reload the index number
	    }

	  }

	  //cout << "hltInfos.size(): " << hltInfos_.size() << endl; 
	
	  //if (hltInfos_.size() > 0 && index != -1)
	  //  cout << "hltNames.size(): " << hltInfos_[index].nHLTPaths() << endl;

	  //cout << "Index: " << index << endl;

	  edm::Handle<edm::TriggerResults> trigResults;
	  try {iEvent.getByLabel(triggerResultsTag_,trigResults);} catch (...) {;}
	  if (trigResults.isValid()) {
	    if (trigResults->wasrun()) nWasRun_++;
	    const bool accept(trigResults->accept());
	    if(verbosity_>0) cout << "   HLT decision: " << accept << endl;
	    rootEvent->setGlobalHLT(accept);
	    if (accept) ++nAccept_;
	    if (trigResults->error() ) nErrors_++;
	  }
	  else
	    {
	      cout << "   HLT results not found!" << endl;;
	      nErrors_++;
	      return;
	    }
	  
	  // decision for each HLT algorithm
	  const unsigned int n(hltNames_.size());
	  std::vector<Bool_t> hltDecision(n, false);
	  for (unsigned int i=0; i!=n; ++i)
	    { 
			  if (trigResults->wasrun(i)) hltInfos_[index].sethltWasRun(i);
	      if (trigResults->error(i) ) hltInfos_[index].sethltErrors(i);
	      if (trigResults->accept(i))
					{
		  			hltInfos_[index].sethltAccept(i);
		  			hltDecision[i]=true;
					}
	    }
	  
	  //cout << "hltDecision.size(): " << hltDecision.size() << endl;
	  rootEvent->setTrigHLT(hltDecision);
	  
	}

	return;
}


void HLTAnalyzer::printStats()
{
	// final printout of accumulated statistics
	
	if(doHLT_)
	{

	  cout << dec << endl;
	  cout << "HLT Summary" << endl;
	  cout << "HLTAnalyzer-Summary " << "---------- Event  Summary ------------\n";
	  cout << "HLTAnalyzer-Summary"
	       << " Events total = " << nEvents_
	       << " wasrun = " << nWasRun_
	       << " passed = " << nAccept_
	       << " errors = " << nErrors_
	       << "\n";
	  
	  cout << endl;
	  cout << "HLTAnalyzer-Summary " << "---------- HLTrig Summary ------------\n";
	  cout << "HLTAnalyzer-Summary "
	       << right << setw(10) << "HLT  Bit#" << " "
	       << right << setw(10) << "WasRun" << " "
	       << right << setw(10) << "Passed" << " "
	       << right << setw(10) << "Errors" << " "
	       << "Name" << "\n";
	  
	  for (unsigned int i = 0; i<hltInfos_.size(); i++) {
	    
	    cat::HLTInfo hltInfo = hltInfos_[i];
	    
	    const unsigned int n(hltInfo.nHLTPaths());
	    
	    for (unsigned int j=0; j!=n; ++j)
	      {
		stringstream s;
		s << hltInfo.RunID();
		cout << "HLTAnalyzer-Summary-RUN "+s.str()
		     << right << setw(10) << j << " "
		     << right << setw(10) << hltInfo.hltWasRun(j) << " "
		     << right << setw(10) << hltInfo.hltAccept(j) << " "
		     << right << setw(10) << hltInfo.hltErrors(j) << " "
		     << hltInfo.hltNames(j) << "\n";
	      }
	  }
	  
	  cout << endl;
	  cout << "HLT Summary end!" << endl;
	  cout << endl;
	}
	
	return;
}

void HLTAnalyzer::copySummary(cat::Run* runInfos)
{
	if(doHLT_)
	{
	  runInfos->setNHLTEvents(nEvents_);
	  runInfos->setNHLTWasRun(nWasRun_);
	  runInfos->setNHLTAccept(nAccept_);
	  runInfos->setNHLTErrors(nErrors_);
	  
	  runInfos->setHLTInputTag(triggerResultsTag_.label()+"_"+triggerResultsTag_.instance()+"_"+triggerResultsTag_.process());
	  cout << triggerResultsTag_.label()+"_"+triggerResultsTag_.instance()+"_"+triggerResultsTag_.process() << endl;
	  
	  // set new hlt container
	  
	  runInfos->setHLTinfos(hltInfos_);
	}

}
