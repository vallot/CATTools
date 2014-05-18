#include "../interface/TrackMETAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

TrackMETAnalyzer::TrackMETAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0), useMC_(false)
{
	metProducer_ = producersNames.getParameter<edm::InputTag>("TrackmetProducer");
	myMETAnalyzer = new RecoMETAnalyzer(producersNames);
}

TrackMETAnalyzer::TrackMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	metProducer_ = producersNames.getParameter<edm::InputTag>("TrackmetProducer");
	useMC_ = myConfig.getUntrackedParameter<bool>("doMETMC");
	myMETAnalyzer = new RecoMETAnalyzer(producersNames,myConfig, verbosity);

}

TrackMETAnalyzer::TrackMETAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	vTrackMETProducer = producersNames.getUntrackedParameter<std::vector<std::string> >("vtrackmetProducer");
	metProducer_ = edm::InputTag(vTrackMETProducer[iter]);
	myMETAnalyzer = new RecoMETAnalyzer(producersNames, myConfig, verbosity);
}

TrackMETAnalyzer::~TrackMETAnalyzer()
{
}

void TrackMETAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootMET)
{

	unsigned int nMETs=0;

	edm::Handle < std::vector <reco::PFMET> > patMETs;
	iEvent.getByLabel(metProducer_, patMETs);
	nMETs = patMETs->size();

	if(verbosity_>1) std::cout << "   Number of MET objects = " << nMETs << "   Label: " << metProducer_.label() << "   Instance: " << metProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nMETs; j++)
	{
		const reco::Candidate* met = 0;	
		met = (const reco::Candidate*) ( & ((*patMETs)[j]) );

		CatMET tempMET = (CatMET) myMETAnalyzer->Process( &( *(met) ) );

		CatTrackMET localMET = CatTrackMET(tempMET);

		localMET.setMETType(2); // 2 = TrackMET

		new( (*rootMET)[j] ) CatTrackMET(localMET);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localMET << endl;
	}
	

}
