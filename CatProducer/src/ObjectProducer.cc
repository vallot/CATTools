#include "../interface/ObjectProducer.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/Handle.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;


ObjectProducer::ObjectProducer(const edm::ParameterSet& iConfig)
{
	myConfig_ = iConfig.getParameter<ParameterSet>("myConfig");
	producersNames_ = iConfig.getParameter<ParameterSet>("producersNames");
}


ObjectProducer::~ObjectProducer()
{
}


// ------------ method called once each job just before starting event loop  ------------
void ObjectProducer::beginJob()
{

	// Load Config parameters	
	verbosity = myConfig_.getUntrackedParameter<int>("verbosity", 0);
	rootFileName_ = myConfig_.getUntrackedParameter<string>("RootFileName","noname.root");
	doHLT = myConfig_.getUntrackedParameter<bool>("doHLT",false);
	doPDFInfo = myConfig_.getUntrackedParameter<bool>("doPDFInfo",false);
	doPrimaryVertex = myConfig_.getUntrackedParameter<bool>("doPrimaryVertex",false);
	runGeneralTracks = myConfig_.getUntrackedParameter<bool>("runGeneralTracks",false);
	doGenJet = myConfig_.getUntrackedParameter<bool>("doGenJet",false);
	doPFJet = myConfig_.getUntrackedParameter<bool>("doPFJet",false);
	doJPTJet = myConfig_.getUntrackedParameter<bool>("doJPTJet",false);
	doMuon = myConfig_.getUntrackedParameter<bool>("doMuon",false);
	doElectron = myConfig_.getUntrackedParameter<bool>("doElectron",false);	
	doPhoton = myConfig_.getUntrackedParameter<bool>("doPhoton",false);	
	doPFMET = myConfig_.getUntrackedParameter<bool>("doPFMET",false);
	doTrackMET = myConfig_.getUntrackedParameter<bool>("doTrackMET",false);
	doTCMET = myConfig_.getUntrackedParameter<bool>("doTCMET",false);
	drawMCTree = myConfig_.getUntrackedParameter<bool>("drawMCTree",false);
	doGenEvent = myConfig_.getUntrackedParameter<bool>("doGenEvent",false);
	doNPGenEvent = myConfig_.getUntrackedParameter<bool>("doNPGenEvent",false);
	doSpinCorrGen = myConfig_.getUntrackedParameter<bool>("doSpinCorrGen",false);
        useEventCounter_ = myConfig_.getUntrackedParameter<bool>("useEventCounter",true);
        filters_ = myConfig_.getUntrackedParameter<std::vector<std::string> >("filters");
	vector<string> defaultVec;
	vGenJetProducer = producersNames_.getUntrackedParameter<vector<string> >("vgenJetProducer",defaultVec);
	vPFJetProducer = producersNames_.getUntrackedParameter<vector<string> >("vpfJetProducer",defaultVec);
	vJPTJetProducer = producersNames_.getUntrackedParameter<vector<string> >("vJPTJetProducer",defaultVec);
	vMuonProducer = producersNames_.getUntrackedParameter<vector<string> >("vmuonProducer",defaultVec);
	vElectronProducer = producersNames_.getUntrackedParameter<vector<string> >("velectronProducer",defaultVec);
	vPhotonProducer = producersNames_.getUntrackedParameter<vector<string> >("vphotonProducer",defaultVec);
 	vPFmetProducer = producersNames_.getUntrackedParameter<vector<string> >("vpfmetProducer",defaultVec);
 	vTrackmetProducer = producersNames_.getUntrackedParameter<vector<string> >("vtrackmetProducer",defaultVec);	

	for(unsigned int s=0;s<vGenJetProducer.size();s++){
		TClonesArray* a;
		vgenJets.push_back(a);
	}

	for(unsigned int s=0;s<vPFJetProducer.size();s++){
		TClonesArray* a;
		vpfJets.push_back(a);
	}

	for(unsigned int s=0;s<vJPTJetProducer.size();s++){
		TClonesArray* a;
		vjptJets.push_back(a);
	}

	for(unsigned int s=0;s<vMuonProducer.size();s++){
		TClonesArray* a;
		vmuons.push_back(a);
	}

	for(unsigned int s=0;s<vElectronProducer.size();s++){
		TClonesArray* a;
		velectrons.push_back(a);
	}
 
        for(unsigned int s=0;s<vPhotonProducer.size();s++){
                TClonesArray* a;
                vphotons.push_back(a);
        }
 
  for(unsigned int s=0; s<vPFmetProducer.size(); s++) {
    TClonesArray* a;
    vPFmets.push_back(a);
  }

 for(unsigned int s=0; s<vTrackmetProducer.size(); s++) {
    TClonesArray* a;
    vTrackmets.push_back(a);
  }


	nTotEvt_ = 0;
	
	// initialize root output file
	rootFile_ = new TFile(rootFileName_.c_str(), "recreate");
	rootFile_->cd();
	if(verbosity>0) cout << "New RootFile " << rootFileName_.c_str() << " is created" << endl;

        tmp_ = new TH1F("EventSummary","EventSummary", filters_.size(),0,filters_.size());

	runInfos_ = new CatRun();
	runTree_ = new TTree("runTree", "Global Run Infos");
	runTree_->Branch ("runInfos", "cat::CatRun", &runInfos_);
	if(verbosity>0) cout << "RunTree is created" << endl;

	rootEvent = 0;
	eventTree_ = new TTree("eventTree", "Event Infos");
	eventTree_->Branch ("Event", "cat::CatEvent", &rootEvent);
	if(verbosity>0) cout << "EventTree is created" << endl;

	if(doHLT)
	{
		if(verbosity>0) cout << "HLT info will be added to rootuple" << endl;
		hltAnalyzer_ = new HLTAnalyzer(producersNames_, myConfig_);
		hltAnalyzer_->setVerbosity(verbosity);
	}

	if(!isRealData_)
	{
		if(verbosity>0) cout << "MC Particles info will be added to rootuple" << endl;
		mcParticles = new TClonesArray("cat::CatMCParticle", 1000);
		eventTree_->Branch ("MCParticles", "TClonesArray", &mcParticles);
	}

	if(!isRealData_ && doGenJet)
	{
		if(verbosity>0) cout << "GenJets info will be added to rootuple (for GenJetStudy)" << endl;
		for(unsigned int s=0; s<vGenJetProducer.size(); s++)
		{
			vgenJets[s] = new TClonesArray("cat::CatGenJet", 1000);
			char name[100];
			sprintf(name,"GenJets_%s",vGenJetProducer[s].c_str());
			eventTree_->Branch (name, "TClonesArray", &vgenJets[s]);
		}
	}
	
	if(doPFJet)
	{
		if(verbosity>0) cout << "PFJets info will be added to rootuple" << endl;
		for(unsigned int s=0;s<vPFJetProducer.size();s++)
		{
			vpfJets[s] = new TClonesArray("cat::CatPFJet", 1000);
			char name[100];
			sprintf(name,"PFJets_%s",vPFJetProducer[s].c_str());
			eventTree_->Branch (name, "TClonesArray", &vpfJets[s]);
		}
	}

	if(doJPTJet)
	{
		if(verbosity>0) cout << "JPT Jets info will be added to rootuple" << endl;
		for(unsigned int s=0;s<vJPTJetProducer.size();s++)
		{
			vjptJets[s] = new TClonesArray("cat::CatJPTJet", 1000);
			char name[100];
			sprintf(name,"JPTJets_%s",vJPTJetProducer[s].c_str());
			eventTree_->Branch (name, "TClonesArray", &vjptJets[s]);
		}
	}
	
	if(doGenEvent)
	{
		if(verbosity>0) cout << "GenEvent info will be added to rootuple" << endl;
		genEvent = new TClonesArray("cat::CatGenEvent", 1000);
		eventTree_->Branch ("GenEvent", "TClonesArray", &genEvent);
	}

	if(doNPGenEvent)
	{
		if(verbosity>0) cout << "NPGenEvent info will be added to rootuple" << endl;
		NPgenEvent = new TClonesArray("cat::CatNPGenEvent", 1000);
		eventTree_->Branch ("NPGenEvent", "TClonesArray", &NPgenEvent);
	}

	if(doSpinCorrGen)
	{
		if(verbosity>0) cout << "SpinCorrelation Gen info will be added to rootuple" << endl;
		spinCorrGen = new TClonesArray("cat::CatSpinCorrGen", 1000);
		eventTree_->Branch ("SpinCorrGen", "TClonesArray", &spinCorrGen);
	}
    
	if(doMuon)
	{
		if(verbosity>0) cout << "Muons info will be added to rootuple" << endl;
		for(unsigned int s=0;s<vMuonProducer.size();s++) {
			vmuons[s] = new TClonesArray("cat::CatMuon", 1000);
			char name[100];
			sprintf(name,"Muons_%s",vMuonProducer[s].c_str());
			eventTree_->Branch (name, "TClonesArray", &vmuons[s]);
		} 
	}
	
	if(doElectron)
	{
		if(verbosity>0) cout << "Electrons info will be added to rootuple" << endl;
		for(unsigned int s=0;s<vElectronProducer.size();s++) {
			velectrons[s] = new TClonesArray("cat::CatElectron", 1000);
			char name[100];
			sprintf(name,"Electrons_%s",vElectronProducer[s].c_str());
			eventTree_->Branch (name, "TClonesArray", &velectrons[s]);
		} 
	}

        if(doPhoton)
        {
                if(verbosity>0) cout << "Photons info will be added to rootuple" << endl;
                for(unsigned int s=0;s<vPhotonProducer.size();s++) {
                        vphotons[s] = new TClonesArray("cat::CatPhoton", 1000);
                        char name[100];
                        sprintf(name,"Photons_%s",vPhotonProducer[s].c_str());
                        eventTree_->Branch (name, "TClonesArray", &vphotons[s]);
                }
        }

	if(doPFMET)
	{
		if(verbosity>0) cout << "ParticleFlowMET info will be added to rootuple" << endl;
    for(unsigned int s=0; s<vPFmetProducer.size(); s++) {
		  vPFmets[s] = new TClonesArray("cat::CatPFMET", 1000);
      char name[100];
			sprintf(name,"PFMET_%s",vPFmetProducer[s].c_str());
  		eventTree_->Branch (name, "TClonesArray", &vPFmets[s]);
	  }
  }


    if(doTrackMET)
        {
                if(verbosity>0) cout << "Track MET info will be added to rootuple" << endl;
    for(unsigned int s=0; s<vTrackmetProducer.size(); s++) {
                  vTrackmets[s] = new TClonesArray("cat::CatTrackMET", 1000);
      char name[100];
                        sprintf(name,"TrackMET_%s",vTrackmetProducer[s].c_str());
                eventTree_->Branch (name, "TClonesArray", &vTrackmets[s]);
          }
  }


	if(doTCMET)
	{
		if(verbosity>0) cout << "Track Corrected MET info will be added to rootuple" << endl;
		TCmet = new TClonesArray("cat::CatMET", 1000);
		eventTree_->Branch ("TCMET", "TClonesArray", &TCmet);
	}
	
	if(doPrimaryVertex)
	{
		if(verbosity>0) cout << "Primary Vertex info will be added to rootuple" << endl;
		primaryVertex = new TClonesArray("cat::CatVertex", 1000);
		eventTree_->Branch ("PrimaryVertex", "TClonesArray", &primaryVertex);
	}

}


// ------------ method called once each job just after ending the event loop  ------------
void ObjectProducer::endJob()
{

	// Trigger Summary Tables
	if(doHLT)
	{	
		cout << "Trigger Summary Tables" << endl;
		hltAnalyzer_->copySummary(runInfos_);
		hltAnalyzer_->printStats();
	}


        runInfos_->setPrePathCounter(tmp_->GetBinContent(1));
        runInfos_->setPostPathCounter(tmp_->GetBinContent(2));
        
	runTree_->Fill();
        //tmp_->Write(); // not saving histogram for consistency, instead we save these numbers in runTree (TJ) 
        std::cout << "Initial number of events: " << tmp_->GetBinContent(1) << std::endl; 
        //Initial number of events could be different from nTotEvt_ if there is any filter in the PAT sequence 
	std::cout << "Total number of events: " << nTotEvt_ << std::endl;
	std::cout << "Closing rootfile " << rootFile_->GetName() << std::endl;
	rootFile_->Write();
	rootFile_->Close();

}

void ObjectProducer::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & setup){
    if(useEventCounter_){
      for(unsigned int i=0;i<filters_.size();++i) {
        std::string name = filters_[i];
        edm::Handle<edm::MergeableCounter> numEventsCounter;
        lumi.getByLabel(name, numEventsCounter);
        if( numEventsCounter.isValid()){
          tmp_->AddBinContent(i+1, numEventsCounter->value);
          tmp_->GetXaxis()->SetBinLabel(i+1,filters_[i].c_str());
        }
      }
    }
}

// ------------ method called to for each event  ------------
void ObjectProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
        isRealData_ = iEvent.isRealData();

	rootFile_->cd();
	nTotEvt_++;
	if( (verbosity>1) || (verbosity>0 && nTotEvt_%10==0 && nTotEvt_<=100)  || (verbosity>0 && nTotEvt_%100==0 && nTotEvt_>100) )
		cout << endl << endl 
			<< "####### ObjectProducer - Cumulated Events " << nTotEvt_
			<< " - Run " << iEvent.id().run() 
			<< " - Event " << iEvent.id().event() 
			<< " #######" << endl;

	// Global Event Infos
	rootEvent = new CatEvent();
	rootEvent->setNb(nTotEvt_);
	rootEvent->setEventId(iEvent.id().event());
	rootEvent->setRunId(iEvent.id().run());

        rootEvent->setLumiBlockId(iEvent.luminosityBlock());

	// do PileUp info

	edm::InputTag PileupSrc_(producersNames_.getParameter<edm::InputTag>("pileUpProducer"));
	Handle<std::vector< PileupSummaryInfo > >  PupInfo;
	iEvent.getByLabel(PileupSrc_, PupInfo);
  
	if (PupInfo.isValid()) {
	  std::vector<PileupSummaryInfo>::const_iterator PVI;
	  
	  // (then, for example, you can do)
	  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	  
	    //std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;

	    rootEvent->setNPu(PVI->getBunchCrossing(),PVI->getPU_NumInteractions());
	    if(PVI->getBunchCrossing() == 0) rootEvent->setNTruePU(PVI->getTrueNumInteractions());
	  }

	  //cout << rootEvent->nPu(0) << endl;
	  //cout << rootEvent->nPu(1) << endl;
	}

	// we need to store some triggerFilter info to be able to emulate triggers on older data
	if (doHLT) {
	  if (verbosity > 1) cout << "should do HLT now..." << endl;
		
	  // get Trigger summary from Event
	  edm::Handle<trigger::TriggerEvent> summary, summary1st, summary2nd, summary3rd, summary4th;
	  edm::InputTag summaryTag1st_("hltTriggerSummaryAOD","",(producersNames_.getParameter < edm::InputTag > ("hltProducer1st")).process());
	  edm::InputTag summaryTag2nd_("hltTriggerSummaryAOD","",(producersNames_.getParameter < edm::InputTag > ("hltProducer2nd")).process());
	  edm::InputTag summaryTag3rd_("hltTriggerSummaryAOD","",(producersNames_.getParameter < edm::InputTag > ("hltProducer3rd")).process());
	  edm::InputTag summaryTag4th_("hltTriggerSummaryAOD","",(producersNames_.getParameter < edm::InputTag > ("hltProducer4th")).process());
	  
	  try { iEvent.getByLabel(summaryTag1st_,summary1st);} catch (...) {;}
	  try { iEvent.getByLabel(summaryTag2nd_,summary2nd);} catch (...) {;}
	  try { iEvent.getByLabel(summaryTag3rd_,summary3rd);} catch (...) {;}
	  try { iEvent.getByLabel(summaryTag4th_,summary4th);} catch (...) {;}

	  //cout << summaryTag1st_ << " " << summaryTag2nd_  << " " << summaryTag3rd_  << " " << summaryTag4th_ << endl;
	  //cout << summaryTag1st_.process() << " " << summaryTag2nd_.process()<< " " << summaryTag3rd_.process()<< " " << summaryTag4th_.process() << endl;
	  //cout << summary1st.isValid() << " " << summary2nd.isValid() << " " << summary3rd.isValid() << " " << summary4th.isValid()<< endl;

	  if (summary1st.isValid()) 
	    summary = summary1st;
	  else if (summary2nd.isValid())
	    summary = summary2nd;
	  else if (summary3rd.isValid())
	    summary = summary3rd;
	  else if (summary4th.isValid())
	    summary = summary4th;
	  else
	    cout << "ObjectProducer::Analyze ERROR: Could not store info for trigger emulation: provided HLTproducerNames are null" << endl;
	    
	  //cout << "summary " << summary << endl;
	  
		if (summary.isValid()) {
	    for (unsigned int i=0; i<summary->sizeFilters(); i++) {
	       if (verbosity > 1) cout << i << " -> " << summary->filterTag(i).label() << endl;
	      
	      // get all trigger objects corresponding to this module.
	      // loop through them and see how many objects match the selection
	      const trigger::Keys& KEYS (summary->filterKeys(i));
	      const int n1(KEYS.size());
	      
	      for (int j=0; j!=n1; ++j) {
					const trigger::TriggerObject& triggerObject( summary-> getObjects().at(KEYS[j]) );
					//cout << "j: " << j << " -> id " << triggerObject.id() << endl;
					//cout << "j: " << j << " -> pt " << triggerObject.pt() << endl;
					//cout << "j: " << j << " -> eta " << triggerObject.eta() << endl;
					//cout << "j: " << j << " -> phi " << triggerObject.phi() << endl;
	      	rootEvent->AddTriggerObject(string(summary->filterTag(i).label()), triggerObject.id(),triggerObject.pt(),triggerObject.eta(),triggerObject.phi());
				}
	    }
	  }
	}
	
	//fastjet density rho
        //FIX ME
	//edm::Handle<double> rho;
	//iEvent.getByLabel("ak5PFJets","rho",rho);
	//rootEvent->setKt6PFJets_rho(*rho);
  
	//density rho for electron isolation (effective area stuff)
  	//edm::Handle<double> rhoIso;
  	//iEvent.getByLabel("kt6PFJetsForIsolation","rho",rhoIso);
  	//rootEvent->setKt6PFJetsForIsolation_rho(*rhoIso);
	
	//if(!isRealData_)
	//{
		//flavorHistory path
		//edm::Handle<unsigned int> flavHist;
		//iEvent.getByLabel("flavorHistoryFilter","",flavHist);
		//rootEvent->setflavorHistoryPath(*flavHist);
	//}
	
	if(runGeneralTracks) // Calculate and fill number of tracks and number of high purity tracks
	{
		// get GeneralTracks collection
		edm::Handle<reco::TrackCollection> tkRef;
		iEvent.getByLabel("generalTracks",tkRef);    
		const reco::TrackCollection* tkColl = tkRef.product();

		if(verbosity>1) std::cout << "Total Number of Tracks " << tkColl->size() << endl;
		rootEvent->setNTracks(tkColl->size());

		int numhighpurity=0;
		reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");

		reco::TrackCollection::const_iterator itk = tkColl->begin();
		reco::TrackCollection::const_iterator itk_e = tkColl->end();
		for(;itk!=itk_e;++itk)
		{
			if(verbosity>1) std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
			if(itk->quality(_trackQuality)) numhighpurity++;
		}

		if(verbosity>1) std::cout << "Total Number of HighPurityTracks " << numhighpurity << endl;
		rootEvent->setNHighPurityTracks(numhighpurity);
	}

	// Trigger
	rootEvent->setGlobalHLT(true);
	if(doHLT)
	{
		if(verbosity>1) std::cout << endl << "Get TriggerResults..." << std::endl;
		//if (nTotEvt_==1) hltAnalyzer_->init(iEvent, rootEvent);
		hltAnalyzer_->process(iEvent, rootEvent);
	}

	// MC Info
	if(!isRealData_)
	{
		if(verbosity>1) cout << endl << "Analysing MC info..." << endl;
		MCAnalyzer* myMCAnalyzer = new MCAnalyzer(myConfig_, producersNames_);
		myMCAnalyzer->SetVerbosity(verbosity);
		if (drawMCTree) myMCAnalyzer->DrawMCTree(iEvent, iSetup, myConfig_, producersNames_);
		if (doPDFInfo ) myMCAnalyzer->PDFInfo(iEvent, rootEvent);
		myMCAnalyzer->ProcessMCParticle(iEvent, mcParticles);
		delete myMCAnalyzer;
	}

	// Get Primary Vertices
	if(doPrimaryVertex)
	{
		if(verbosity>1) cout << endl << "Analysing primary vertices collection..." << endl;
		VertexAnalyzer* myVertexAnalyzer = new VertexAnalyzer(producersNames_, verbosity);
		myVertexAnalyzer->Process(iEvent, primaryVertex);
		delete myVertexAnalyzer;
	}

	// GenJet
	if(!isRealData_ && doGenJet)
	{
		if(verbosity>1) cout << endl << "Analysing GenJets collection ..." << endl;
		for(unsigned int s=0; s<vGenJetProducer.size(); s++)
		{
			GenJetAnalyzer* myGenJetAnalyzer = new GenJetAnalyzer(producersNames_, s, myConfig_, verbosity);
			myGenJetAnalyzer->Process(iEvent, vgenJets[s]);
			delete myGenJetAnalyzer;
		}
	}

	// PFJet
	if(doPFJet)
	{
		if(verbosity>1) cout << endl << "Analysing PFjets collection..." << endl;
		for(unsigned int s=0;s<vPFJetProducer.size();s++){
			PFJetAnalyzer* myPFJetAnalyzer = new PFJetAnalyzer(producersNames_, s,  myConfig_, verbosity);
			myPFJetAnalyzer->Process(iEvent, vpfJets[s], iSetup);
			delete myPFJetAnalyzer;
		}
	}

	// JPT Jets
	if(doJPTJet)
	{
		if(verbosity>1) cout << endl << "Analysing JPT jets collection..." << endl;
		for(unsigned int s=0;s<vJPTJetProducer.size();s++){
			JPTJetAnalyzer* myJPTJetAnalyzer = new JPTJetAnalyzer(producersNames_, s,  myConfig_, verbosity);
			myJPTJetAnalyzer->Process(iEvent, vjptJets[s], iSetup);
			delete myJPTJetAnalyzer;
		}
	}

	// GenEvent
	if(doGenEvent)
	{
		if(verbosity>1) cout << endl << "Analysing GenEvent collection..." << endl;
		GenEventAnalyzer* myGenEventAnalyzer = new GenEventAnalyzer(producersNames_, myConfig_, verbosity);
		myGenEventAnalyzer->Process(iEvent, genEvent);
		delete myGenEventAnalyzer;
	}

	// NPGenEvent
	if(doNPGenEvent)
	{
		if(verbosity>1) cout << endl << "Analysing NPGenEvent collection..." << endl;
		NPGenEventAnalyzer* myNPGenEventAnalyzer = new NPGenEventAnalyzer(producersNames_, myConfig_, verbosity);
		if(verbosity>1) cout << endl << "Analysing NPGenEvent collection..." << endl;
		myNPGenEventAnalyzer->Process(iEvent, NPgenEvent);
		if(verbosity>1) cout << endl << "Analysing NPGenEvent collection..." << endl;
		delete myNPGenEventAnalyzer;
	}

	// SpinCorrelation Gen
	if(doSpinCorrGen)
	{
		if(verbosity>1) cout << endl << "Analysing SpinCorrGen collection..." << endl;
		SpinCorrGenAnalyzer* mySpinCorrGenAnalyzer = new SpinCorrGenAnalyzer(producersNames_, myConfig_, verbosity);
		mySpinCorrGenAnalyzer->Process(iEvent, spinCorrGen);
		delete mySpinCorrGenAnalyzer;
	}

	// Muons
	if(doMuon)
	{
		if(verbosity>1) cout << endl << "Analysing muon collection..." << endl;
		for(unsigned int s=0;s<vMuonProducer.size();s++){
		  MuonAnalyzer* myMuonAnalyzer = new MuonAnalyzer(producersNames_, s, myConfig_, verbosity);
		  myMuonAnalyzer->Process(iEvent, vmuons[s]);
		  delete myMuonAnalyzer;
		}
	}	

	// Electrons
	if(doElectron)
	{
		if(verbosity>1) cout << endl << "Analysing electrons collection..." << endl;
		for(unsigned int s=0;s<vElectronProducer.size();s++){
		  ElectronAnalyzer* myElectronAnalyzer = new ElectronAnalyzer(producersNames_, s, myConfig_, verbosity);
		  myElectronAnalyzer->Process(iEvent, velectrons[s], iSetup);
		  delete myElectronAnalyzer;
		}
	}	

        // Photons
        if(doPhoton)
        {
                if(verbosity>1) cout << endl << "Analysing photons collection..." << endl;
                for(unsigned int s=0;s<vPhotonProducer.size();s++){
                  PhotonAnalyzer* myPhotonAnalyzer = new PhotonAnalyzer(producersNames_, s, myConfig_, verbosity);
                  myPhotonAnalyzer->Process(iEvent, vphotons[s], iSetup);
                  delete myPhotonAnalyzer;
                }
        }

	if(doPFMET)
	{
		if(verbosity>1) cout << endl << "Analysing ParticleFlow Missing Et..." << endl;
    for(unsigned int s=0; s<vPFmetProducer.size(); s++) {
      PFMETAnalyzer* myPFMETAnalyzer = new PFMETAnalyzer(producersNames_, s, myConfig_, verbosity);
      myPFMETAnalyzer->Process(iEvent, vPFmets[s]);
      delete myPFMETAnalyzer;
    }
	}

 if(doTrackMET)
        {
                if(verbosity>1) cout << endl << "Analysing Track Missing Et..." << endl;
    for(unsigned int s=0; s<vTrackmetProducer.size(); s++) {
      TrackMETAnalyzer* myTrackMETAnalyzer = new TrackMETAnalyzer(producersNames_, s, myConfig_, verbosity);
      myTrackMETAnalyzer->Process(iEvent, vTrackmets[s]);
      delete myTrackMETAnalyzer;
    }
        }



	if(doTCMET)
	{
		if(verbosity>1) cout << endl << "Analysing Track Corrected Missing Et..." << endl;
		TCMETAnalyzer* myMETAnalyzer = new TCMETAnalyzer(producersNames_, myConfig_, verbosity);
		myMETAnalyzer->Process(iEvent, TCmet);
		delete myMETAnalyzer;
	}
	
	// Associate recoParticles to mcParticles
	if(!isRealData_)
	{
		MCAssociator* myMCAssociator = new MCAssociator(producersNames_, verbosity);
		myMCAssociator->init(iEvent, mcParticles);
		if(doPFJet && vpfJets.size() > 0) myMCAssociator->process(vpfJets[0]);
		if(doMuon && vmuons.size() > 0) myMCAssociator->process(vmuons[0]);
		if(doElectron && velectrons.size() > 0) myMCAssociator->process(velectrons[0]);
		//if(verbosity>2 && doMuon && vmuons.size() > 0) myMCAssociator->printParticleAssociation(vmuons[0]);
		//if(verbosity>2 && doElectron && velectrons.size() > 0) myMCAssociator->printParticleAssociation(velectrons[0]);
		//if(verbosity>2 && doPhoton) myMCAssociator->printParticleAssociation(photons);
		delete myMCAssociator;
	}


	if(verbosity>1) cout << endl << "Filling rootuple..." << endl;
	eventTree_->Fill();
	if(verbosity>1) cout << endl << "Deleting objects..." << endl;
	delete rootEvent;
	if(!isRealData_) (*mcParticles).Delete();
	if(!isRealData_ && doGenJet)
	{
		for(unsigned int s=0;s<vGenJetProducer.size();s++)
		{
			(*vgenJets[s]).Delete();
		}		
	}
	if(doPFJet){
		for(unsigned int s=0;s<vPFJetProducer.size();s++){
			(*vpfJets[s]).Delete();
		}
	}
	if(doJPTJet){
		for(unsigned int s=0;s<vJPTJetProducer.size();s++){
			(*vjptJets[s]).Delete();
		}
	}
	if(doMuon){
		for(unsigned int s=0;s<vMuonProducer.size();s++){
			(*vmuons[s]).Delete();
		}
	}
	if(doElectron){
		for(unsigned int s=0;s<vElectronProducer.size();s++){
			(*velectrons[s]).Delete();
		}
	}
        if(doPhoton){
                for(unsigned int s=0;s<vPhotonProducer.size();s++){
                        (*vphotons[s]).Delete();
                }
        }
	if(doPFMET) {
    for(unsigned int s=0; s<vPFmetProducer.size(); s++) {
      (*vPFmets[s]).Delete();
    }
  }
 if(doTrackMET) {
    for(unsigned int s=0; s<vTrackmetProducer.size(); s++) {
      (*vTrackmets[s]).Delete();
    }
  }

	if(doTCMET) (*TCmet).Delete();
	if(doGenEvent) (*genEvent).Delete();
	if(doNPGenEvent) (*NPgenEvent).Delete();
	if(doSpinCorrGen) (*spinCorrGen).Delete();
	if(doPrimaryVertex) (*primaryVertex).Delete();
	if(verbosity>0) cout << endl;
}
