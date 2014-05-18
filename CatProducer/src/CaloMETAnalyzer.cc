#include "../interface/CaloMETAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

CaloMETAnalyzer::CaloMETAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0), useMC_(false)
{
	metProducer_ = producersNames.getParameter<edm::InputTag>("CalometProducer");
	myMETAnalyzer = new METAnalyzer(producersNames);
}

CaloMETAnalyzer::CaloMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	metProducer_ = producersNames.getParameter<edm::InputTag>("CalometProducer");
	useMC_ = myConfig.getUntrackedParameter<bool>("doMETMC");
	myMETAnalyzer = new METAnalyzer(producersNames,myConfig, verbosity);

}

CaloMETAnalyzer::~CaloMETAnalyzer()
{
}

void CaloMETAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootMET)
{

	unsigned int nMETs=0;

	edm::Handle < std::vector <reco::CaloMET> > recoMETs;
	edm::Handle < std::vector <pat::MET> > patMETs;
	iEvent.getByLabel(metProducer_, patMETs);
	nMETs = patMETs->size();
	
	if(verbosity_>1) std::cout << "   Number of MET objects = " << nMETs << "   Label: " << metProducer_.label() << "   Instance: " << metProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nMETs; j++)
	{
		const reco::Candidate* met = 0;	
		met = (const reco::Candidate*) ( & ((*patMETs)[j]) );

		CatMET tempMET = (CatMET) myMETAnalyzer->Process( &( *(met) ) );

		CatCaloMET localMET = CatCaloMET(tempMET);

		localMET.setMETType(1); // 1 = CaloMET

		const pat::MET *patMET = dynamic_cast<const pat::MET*>(&*met);

			localMET.setCaloMETFraction(
				patMET->maxEtInEmTowers()
				,patMET->maxEtInHadTowers()
				,patMET->hadEtInHO()
				,patMET->hadEtInHB()
				,patMET->hadEtInHF()
				,patMET->hadEtInHE()
				,patMET->emEtInEB()
				,patMET->emEtInEE()
				,patMET->emEtInHF()
				,patMET->etFractionHadronic()
				,patMET->emEtFraction()
				,patMET->metSignificance()
				,patMET->CaloMETInpHF()
				,patMET->CaloMETInmHF()
				,patMET->CaloSETInpHF()
				,patMET->CaloSETInmHF()
				,patMET->CaloMETPhiInpHF()
				,patMET->CaloMETPhiInmHF()
				);
			
		new( (*rootMET)[j] ) CatCaloMET(localMET);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localMET << endl;
	}

}
