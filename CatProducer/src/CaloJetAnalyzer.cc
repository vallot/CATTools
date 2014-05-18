#include "../interface/CaloJetAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

CaloJetAnalyzer::CaloJetAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	caloJetProducer_ = producersNames.getParameter<edm::InputTag>("caloJetProducer");
	myJetAnalyzer = new JetAnalyzer();
}

CaloJetAnalyzer::CaloJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	caloJetProducer_ = producersNames.getParameter<edm::InputTag>("caloJetProducer");
	myJetAnalyzer = new JetAnalyzer(verbosity);
}

CaloJetAnalyzer::CaloJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	caloJetProducer_ = producersNames.getParameter<edm::InputTag>("caloJetProducer");
	doCaloJetId_ = myConfig.getUntrackedParameter<bool>("doCaloJetId");
	myJetAnalyzer = new JetAnalyzer(myConfig, verbosity);
}

CaloJetAnalyzer::CaloJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	vCaloJetProducer = producersNames.getUntrackedParameter<std::vector<std::string> >("vcaloJetProducer");
	caloJetProducer_ = edm::InputTag(vCaloJetProducer[iter]);
	doCaloJetId_ = myConfig.getUntrackedParameter<bool>("doCaloJetId");
	myJetAnalyzer = new JetAnalyzer(myConfig, verbosity);
}

CaloJetAnalyzer::~CaloJetAnalyzer()
{
}

void CaloJetAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootJets, const edm::EventSetup& iSetup)
{

	unsigned int nJets=0;

	// check if the jet is of the good type
	std::string jetType = "BASIC";
	if( caloJetProducer_.label()=="kt4CaloJets"
		|| caloJetProducer_.label()=="kt6CaloJets"
		|| caloJetProducer_.label()=="iterativeCone5CaloJets"
		|| caloJetProducer_.label()=="sisCone5CaloJets"
		|| caloJetProducer_.label()=="sisCone7CaloJets"
      || caloJetProducer_.label()=="ak5CaloJets"
      || caloJetProducer_.label()=="ak7CaloJets"
	) jetType="CALO";

	edm::Handle < std::vector <pat::Jet> > patJets;
	iEvent.getByLabel(caloJetProducer_, patJets);
	nJets = patJets->size();
		
	if(verbosity_>1) std::cout << "   Number of jets = " << nJets << "   Label: " << caloJetProducer_.label() << "   Instance: " << caloJetProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nJets; j++)
	{
		const reco::Jet* jet = 0;	
		jet = (const reco::Jet*) ( & ((*patJets)[j]) );
		if( (*patJets)[j].isCaloJet() ) jetType="CALO";
		
		// Call JetAnalyzer to fill the basic Jet Properties
		CatJet tempJet = (CatJet) myJetAnalyzer->Process( &( *(jet) ), iSetup);

		CatCaloJet localJet = CatCaloJet(tempJet);

		localJet.setJetType(1); // 1 = CaloJet
		localJet.setN90(jet->nCarrying(0.9));
		localJet.setN60(jet->nCarrying(0.6));
		localJet.setetaetaMoment(jet->etaetaMoment());
		localJet.setphiphiMoment(jet->phiphiMoment());

		// Some specific methods to pat::Jet
		const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);
			
		localJet.setChargedMultiplicity(patJet->associatedTracks().size()) ;
		localJet.setEcalEnergyFraction(patJet->emEnergyFraction());
		localJet.setHcalEnergyFraction(patJet->energyFractionHadronic());
		localJet.setMaxEInEmTowers(patJet->maxEInEmTowers());
		localJet.setMaxEInHadTowers(patJet->maxEInHadTowers());
		localJet.setTowersArea(patJet->towersArea());

		if(doCaloJetId_)
		{
			localJet.setfHPD(patJet->jetID().fHPD);
			localJet.setfRBX(patJet->jetID().fRBX);
			localJet.setn90Hits(patJet->jetID().n90Hits);
			localJet.setnHCALTowers(patJet->jetID().nHCALTowers);
			localJet.setnECALTowers(patJet->jetID().nECALTowers);
		} //end of if(doCaloJetId_)

		new( (*rootJets)[j] ) CatCaloJet(localJet);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localJet << endl;
		
	}

}
