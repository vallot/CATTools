#include "../interface/TCMETAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

TCMETAnalyzer::TCMETAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0), useMC_(false)
{
	metProducer_ = producersNames.getParameter<edm::InputTag>("TCmetProducer");
	myMETAnalyzer = new METAnalyzer(producersNames);
}

TCMETAnalyzer::TCMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	metProducer_ = producersNames.getParameter<edm::InputTag>("TCmetProducer");
	useMC_ = myConfig.getUntrackedParameter<bool>("doMETMC");
	myMETAnalyzer = new METAnalyzer(producersNames,myConfig, verbosity);

}

TCMETAnalyzer::~TCMETAnalyzer()
{
}

void TCMETAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootMET)
{

	unsigned int nMETs=0;

	edm::Handle < std::vector <pat::MET> > patMETs;
	iEvent.getByLabel(metProducer_, patMETs);
	nMETs = patMETs->size();

	if(verbosity_>1) std::cout << "   Number of MET objects = " << nMETs << "   Label: " << metProducer_.label() << "   Instance: " << metProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nMETs; j++)
	{
		const reco::Candidate* met = 0;	
		met = (const reco::Candidate*) ( & ((*patMETs)[j]) );

		CatMET tempMET = (CatMET) myMETAnalyzer->Process( &( *(met) ) );

		CatMET localMET = CatMET(tempMET);

		localMET.setMETType(3); // 3 = TCMET
				
		new( (*rootMET)[j] ) CatMET(localMET);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localMET << endl;
	}

}
