#include "../interface/GenJetAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

GenJetAnalyzer::GenJetAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	genJetProducer_ = producersNames.getParameter<edm::InputTag>("genJetProducer");
}

GenJetAnalyzer::GenJetAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	genJetProducer_ = producersNames.getParameter<edm::InputTag>("genJetProducer");
}

GenJetAnalyzer::GenJetAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	genJetProducer_ = producersNames.getParameter<edm::InputTag>("genJetProducer");
}

GenJetAnalyzer::GenJetAnalyzer(const edm::ParameterSet& producersNames, int iter, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	vGenJetProducer = producersNames.getUntrackedParameter<std::vector<std::string> >("vgenJetProducer");
	genJetProducer_ = edm::InputTag(vGenJetProducer[iter]);
}

GenJetAnalyzer::~GenJetAnalyzer()
{
}

void GenJetAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootGenJets)
{
	// check if the genJet is of the good type
	std::string jetType = "BASIC";
	if( genJetProducer_.label()=="kt4GenJets" || genJetProducer_.label()=="kt6GenJets" || genJetProducer_.label()=="ak5GenJets"
	|| genJetProducer_.label()=="ak7GenJets" || genJetProducer_.label()=="ak5GenJetsNoE" || genJetProducer_.label()=="ak5GenJetsNoNu" || genJetProducer_.label()=="ak5GenJetsNoMuNoNu")
		jetType="CALO";

	edm::Handle < edm::View <reco::GenJet> > recoGenJets;
	iEvent.getByLabel(genJetProducer_, recoGenJets);
	
	unsigned int nJets = recoGenJets->size();

	if(verbosity_>1) std::cout << "   Number of jets = " << nJets << "   Label: " << genJetProducer_.label() << "   Instance: " << genJetProducer_.instance() << std::endl;

	for (unsigned int j=0; j<nJets; j++)
	{
		const reco::GenJet* genJet = 0;	
		if( jetType=="CALO" ) genJet = (const reco::GenJet*) ( & ((*recoGenJets)[j]) );
			
		// Call JetAnalyzer to fill the basic Jet Properties
//		cat::CatJet tempJet = (cat::CatJet) myJetAnalyzer->Process( &( *(genJet) ));

		cat::CatGenJet localGenJet(
			genJet->px()
			,genJet->py()
			,genJet->pz()
			,genJet->energy()
			,genJet->vx()
			,genJet->vy()
			,genJet->vz()
			,genJet->pdgId()
			,genJet->charge()
		);

		localGenJet.setNConstituents(genJet->nConstituents());
		localGenJet.setMaxDistance(genJet->maxDistance());
		localGenJet.setN90(genJet->nCarrying(0.9));
		localGenJet.setN60(genJet->nCarrying(0.6));
		localGenJet.setetaetaMoment(genJet->etaetaMoment());
		localGenJet.setphiphiMoment(genJet->phiphiMoment());
		localGenJet.setEMEnergy(genJet->emEnergy());
		localGenJet.setHadEnergy(genJet->hadEnergy());
		localGenJet.setInvisibleEnergy(genJet->invisibleEnergy());

                bool isBHadron = false;
                bool isCHadron = false;
                cat::CatMCParticle BHad;
                cat::CatMCParticle CHad;

                std::vector <const reco::GenParticle*> mcparts = genJet->getGenConstituents();

                for (unsigned i = 0; i < mcparts.size (); i++) {
                  const GenParticle* mcpart = mcparts[i];
                  const reco::Candidate* lastB = lastBHadron(*mcpart);
                  if( lastB ) {
                    isBHadron = true;
                    cat::CatMCParticle tmp( lastB->px(), lastB->py(), lastB->pz(), lastB->energy() );
                    BHad = tmp;
                    break;
                  }
                }

                for (unsigned i = 0; i < mcparts.size (); i++) {
                  if( isBHadron ) break; //no need to loop over again, this is b-jet!
                  const GenParticle* mcpart = mcparts[i];
                  const reco::Candidate* lastC = lastCHadron(*mcpart);
                  if( lastC ) {
                    isCHadron = true;
                    cat::CatMCParticle tmp( lastC->px(), lastC->py(), lastC->pz(), lastC->energy() );
                    CHad = tmp;
                    break;
                  }
                }
                 
                if( isBHadron ) localGenJet.setBHadron(BHad); //if B-Hadron matched, always assign B-Hadron
                if( isCHadron ) localGenJet.setCHadron(CHad); //if only no B-Hadron matched, assign C-Hadron
                
				
		new( (*rootGenJets)[j] ) cat::CatGenJet(localGenJet);
		if(verbosity_>2) cout << "   ["<< setw(3) << j << "] " << localGenJet << endl;
	}
}

std::vector<const reco::Candidate *> GenJetAnalyzer::getAncestors(const reco::Candidate &c)
{
  vector<const reco::Candidate *> moms;
  if( c.numberOfMothers() == 1 ) {
    const Candidate * dau = &c;
    const Candidate * mom = c.mother();
    while ( dau->numberOfMothers() == 1) {
      moms.push_back( dau );
      dau = mom ;
      mom = dau->mother();
    }
  }
  return moms;
}

bool GenJetAnalyzer::hasBottom(const reco::Candidate &c)
{
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;
}

bool GenJetAnalyzer::hasCharm(const reco::Candidate &c)
{
  int code1;
  int code2;
  bool tmpHasCharm = false;
  code1 = (int)( ( abs(c.pdgId() ) / 100)%10 );
  code2 = (int)( ( abs(c.pdgId() ) /1000)%10 );
  if ( code1 == 4 || code2 == 4) tmpHasCharm = true;
  return tmpHasCharm;
}

bool GenJetAnalyzer::decayFromBHadron(const Candidate & c)
{
   bool isFromB = false;
   vector<const Candidate *> allParents = getAncestors( c );
   for( vector<const Candidate *>::const_iterator aParent = allParents.begin();
                                                  aParent != allParents.end();
                                                  aParent ++ )
     {
         if( hasBottom(**aParent) ) isFromB = true;
/*
 cout << " particle Parent is " << (*aParent)->status()
 << " type " << (*aParent)->pdgId()
 << " pt= " << (*aParent)->pt()
 << " isB = " << isFromB
 << endl;
*/
     }
   return isFromB;
}

bool GenJetAnalyzer::decayFromCHadron(const Candidate & c)
{
  bool isFromC = false;
  vector<const Candidate *> allParents = getAncestors( c );
  for( vector<const Candidate *>::const_iterator aParent = allParents.begin();
                                                 aParent != allParents.end();
                                                 aParent ++ )
  {
    if( hasCharm(**aParent) ) isFromC = true;
/*
cout << " particle Parent is " << (*aParent)->status()
<< " type " << (*aParent)->pdgId()
<< " pt=" << (*aParent)->pt()
<< " isC = " << isFromC
<< endl;
*/
   }
   return isFromC;
}

const Candidate* GenJetAnalyzer::lastBHadron(const Candidate & c)
{
   const Candidate * out = 0;
   
   vector<const Candidate *> allParents = getAncestors( c );
   for( vector<const Candidate *>::const_iterator aParent = allParents.begin();
                                                  aParent != allParents.end();
                                                  aParent ++ )
     {
         if( hasBottom(**aParent) ) out = *aParent;
         
     }
   return out;
}

const Candidate* GenJetAnalyzer::lastCHadron(const Candidate & c)
{
   const Candidate * out = 0;

   vector<const Candidate *> allParents = getAncestors( c );
   for( vector<const Candidate *>::const_iterator aParent = allParents.begin();
                                                  aParent != allParents.end();
                                                  aParent ++ )
     {
         if( hasCharm(**aParent) ) out = *aParent;

     }
   return out;
}
