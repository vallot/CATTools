#include "../interface/JPTJetAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

JPTJetAnalyzer::JPTJetAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	JPTJetProducer_ = producersNames.getParameter<edm::InputTag>("JPTJetProducer");
	myJetAnalyzer = new JetAnalyzer();
}

JPTJetAnalyzer::JPTJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	JPTJetProducer_ = producersNames.getParameter<edm::InputTag>("JPTJetProducer");
	myJetAnalyzer = new JetAnalyzer(verbosity);
}

JPTJetAnalyzer::JPTJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	JPTJetProducer_ = producersNames.getParameter<edm::InputTag>("JPTJetProducer");
	doJPTJetId_ = myConfig.getUntrackedParameter<bool>("doJPTJetId");
	myJetAnalyzer = new JetAnalyzer(myConfig, verbosity);
}

JPTJetAnalyzer::JPTJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	vJPTJetProducer = producersNames.getUntrackedParameter<std::vector<std::string> >("vJPTJetProducer");
	JPTJetProducer_ = edm::InputTag(vJPTJetProducer[iter]);
	doJPTJetId_ = myConfig.getUntrackedParameter<bool>("doJPTJetId");
	myJetAnalyzer = new JetAnalyzer(myConfig, verbosity);
}

JPTJetAnalyzer::~JPTJetAnalyzer()
{
}

void JPTJetAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootJets, const edm::EventSetup& iSetup)
{

	unsigned int nJets=0;

	// check if the jet is of the good type
	std::string jetType = "BASIC";
	
	edm::Handle < std::vector <pat::Jet> > patJets;
	iEvent.getByLabel(JPTJetProducer_, patJets);
	nJets = patJets->size();
		
	if(verbosity_>1) std::cout << "   Number of jets = " << nJets << "   Label: " << JPTJetProducer_.label() << "   Instance: " << JPTJetProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nJets; j++)
	{
		const reco::Jet* jet = 0;	
	
		jet = (const reco::Jet*) ( & ((*patJets)[j]) );
		if( (*patJets)[j].isJPTJet() ) jetType="JPT";
			
		// Call JetAnalyzer to fill the basic Jet Properties
		CatJet tempJet = (CatJet) myJetAnalyzer->Process( &( *(jet) ), iSetup);

		CatJPTJet localJet = CatJPTJet(tempJet);

		localJet.setJetType(3); // 3 = JPTJet
		localJet.setN90(jet->nCarrying(0.9));
		localJet.setN60(jet->nCarrying(0.6));
		localJet.setetaetaMoment(jet->etaetaMoment());
		localJet.setphiphiMoment(jet->phiphiMoment());

		// Some specific methods to pat::Jet
		const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);
			
		localJet.setEcalEnergyFraction(patJet->emEnergyFraction());
		localJet.setHcalEnergyFraction(patJet->energyFractionHadronic());
		localJet.setMaxEInEmTowers(patJet->maxEInEmTowers());
		localJet.setMaxEInHadTowers(patJet->maxEInHadTowers());
		localJet.setTowersArea(patJet->towersArea());
		localJet.setChargedMultiplicity(patJet->chargedMultiplicity()) ;
		localJet.setchargedHadronEnergyFraction(patJet->chargedHadronEnergyFraction());
		
		if(doJPTJetId_)
		{
			localJet.setfHPD(patJet->jetID().fHPD);
			localJet.setfRBX(patJet->jetID().fRBX);
			localJet.setn90Hits(patJet->jetID().n90Hits);
			localJet.setnHCALTowers(patJet->jetID().nHCALTowers);
			localJet.setnECALTowers(patJet->jetID().nECALTowers);
		} //end of if(doJPTJetId_)	
			
		new( (*rootJets)[j] ) CatJPTJet(localJet);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localJet << endl;
		
	}
}
