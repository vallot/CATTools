#include "../interface/PFJetAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

PFJetAnalyzer::PFJetAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	pfJetProducer_ = producersNames.getParameter<edm::InputTag>("pfJetProducer");
	myJetAnalyzer = new JetAnalyzer();
}

PFJetAnalyzer::PFJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	pfJetProducer_ = producersNames.getParameter<edm::InputTag>("pfJetProducer");
	myJetAnalyzer = new JetAnalyzer(verbosity);
}

PFJetAnalyzer::PFJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	pfJetProducer_ = producersNames.getParameter<edm::InputTag>("pfJetProducer");
	myJetAnalyzer = new JetAnalyzer(myConfig, verbosity);
}

PFJetAnalyzer::PFJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	vPFJetProducer = producersNames.getUntrackedParameter<std::vector<std::string> >("vpfJetProducer");
	pfJetProducer_ = edm::InputTag(vPFJetProducer[iter]);
	myJetAnalyzer = new JetAnalyzer(myConfig, verbosity);
}

PFJetAnalyzer::~PFJetAnalyzer()
{
}

void PFJetAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootJets, const edm::EventSetup& iSetup)
{

	unsigned int nJets=0;

	// check if the jet is of the good type
	std::string jetType = "BASIC";
	if( pfJetProducer_.label()=="kt4PFJets"
		|| pfJetProducer_.label()=="kt6PFJets"
		|| pfJetProducer_.label()=="iterativeCone5PFJets"
		|| pfJetProducer_.label()=="sisCone5PFJets"
		|| pfJetProducer_.label()=="sisCone7PFJets"
		|| pfJetProducer_.label()=="ak5PFJets"
		|| pfJetProducer_.label()=="ak7PFJets"
	) jetType="PF";

	edm::Handle < std::vector <pat::Jet> > patJets;
	iEvent.getByLabel(pfJetProducer_, patJets);
	nJets = patJets->size();
	
		
	if(verbosity_>1) std::cout << "   Number of jets = " << nJets << "   Label: " << pfJetProducer_.label() << "   Instance: " << pfJetProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nJets; j++)
	{
		const reco::Jet* jet = 0;	
		jet = (const reco::Jet*) ( & ((*patJets)[j]) );
		if( (*patJets)[j].isPFJet() ) jetType="PF";
			
		// Call JetAnalyzer to fill the basic Jet Properties
		cat::Jet tempJet = myJetAnalyzer->Process( &( *(jet) ), iSetup);
		
		cat::PFJet localJet = cat::PFJet(tempJet);

		localJet.setJetType(2); // 2 = PFJet

		// Some specific methods to pat::Jet
		const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);
			
		localJet.setChargedHadronEnergyFraction(patJet->chargedHadronEnergyFraction());
		localJet.setNeutralHadronEnergyFraction(patJet->neutralHadronEnergyFraction());
		localJet.setChargedEmEnergyFraction(patJet->chargedEmEnergyFraction());
		localJet.setChargedMuEnergyFraction(patJet->chargedMuEnergyFraction());
		localJet.setNeutralEmEnergyFraction(patJet->neutralEmEnergyFraction());		
		localJet.setHFHadronEnergyFraction(patJet->HFHadronEnergyFraction());
		localJet.setHFEMEnergyFraction(patJet->HFEMEnergyFraction());		
		localJet.setChargedMultiplicity(patJet->chargedMultiplicity());
		localJet.setNeutralMultiplicity(patJet->neutralMultiplicity());
		localJet.setMuonMultiplicity(patJet->muonMultiplicity());		
		localJet.setHFHadronMultiplicity(patJet->HFHadronMultiplicity());
		localJet.setHFEMMultiplicity(patJet->HFEMMultiplicity());

		new( (*rootJets)[j] ) cat::PFJet(localJet);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localJet << endl;
		
	}

}
