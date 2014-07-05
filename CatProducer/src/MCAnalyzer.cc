#include "../interface/MCAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;


MCAnalyzer::MCAnalyzer():
	verbosity_(0)
	,doElectronMC_(false)
	,electronMC_etaMax_(99999.)
	,electronMC_ptMin_(99999.)
	,doMuonMC_(false)
	,muonMC_etaMax_(99999.)
	,muonMC_ptMin_(99999.)
	,doUnstablePartsMC_(false)
	,signalGenerator_("noname")
{
}



MCAnalyzer::MCAnalyzer(const edm::ParameterSet& config, const edm::ParameterSet& producersNames) : verbosity_(0)
{
	
	doElectronMC_ = config.getUntrackedParameter<bool>("doElectronMC", false);
	electronMC_etaMax_ = config.getParameter<double>("electronMC_etaMax");
	electronMC_ptMin_ = config.getParameter<double>("electronMC_ptMin");
	
	doMuonMC_ = config.getUntrackedParameter<bool>("doMuonMC", false);
	muonMC_etaMax_ = config.getParameter<double>("muonMC_etaMax");
	muonMC_ptMin_ = config.getParameter<double>("muonMC_ptMin");

	doJetMC_ = config.getUntrackedParameter<bool>("doJetMC", false);
	jetMC_etaMax_ = config.getParameter<double>("jetMC_etaMax");
	jetMC_ptMin_ = config.getParameter<double>("jetMC_ptMin");

	doMETMC_ = config.getUntrackedParameter<bool>("doMETMC", false);

	doUnstablePartsMC_ = config.getUntrackedParameter<bool>("doUnstablePartsMC", false);

	signalGenerator_ = config.getUntrackedParameter<string>("signalGenerator","noname");
	genParticlesProducer_ = producersNames.getParameter<edm::InputTag>("genParticlesProducer");
}



MCAnalyzer::~MCAnalyzer()
{
}



void MCAnalyzer::DrawMCTree(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::ParameterSet& config, const edm::ParameterSet& producersNames)
{
	cout << endl << " ----- ParticleTreeDrawer -----" << endl;
	ParticleTreeDrawer* ptd = new ParticleTreeDrawer(config, producersNames);
	ptd->analyze( iEvent, iSetup );
	delete ptd;
}



void MCAnalyzer::PDFInfo(const edm::Event& iEvent, cat::Event* rootEvent)
{
	if(verbosity_>1) cout << endl << "   Process PDF Infos..." << endl;
	edm::Handle<GenEventInfoProduct> genEvtInfo;
	iEvent.getByLabel( "generator", genEvtInfo );
	typedef gen::PdfInfo PDF;
	const PDF *pdfInfo = genEvtInfo->pdf();
	if (genEvtInfo->hasPDF() && verbosity_>1)
	{
		cout << "   First incoming parton:  flavour=" << pdfInfo->id.first << " x1 = " << pdfInfo->x.first << endl;
		cout << "   Second incoming parton: flavour=" << pdfInfo->id.second << " x2 = " << pdfInfo->x.second << endl;
		cout << "   Factorization Scale Q = " << pdfInfo->scalePDF << endl;
	}
	rootEvent->setIdParton1(pdfInfo->id.first);
	rootEvent->setXParton1(pdfInfo->x.first);		
	rootEvent->setIdParton2(pdfInfo->id.second);
	rootEvent->setXParton2(pdfInfo->x.second);
	rootEvent->setFactorizationScale(pdfInfo->scalePDF);	
	
}



void MCAnalyzer::ProcessMCParticle(const edm::Event& iEvent, TClonesArray* rootMCParticles)
{
	// Fill TCloneArrays with preselected MC Electrons, Muons  and with the primary decaying particles
	if(verbosity_>1) cout << endl << "   Process MC Particles..." << endl;
	edm::Handle <reco::GenParticleCollection> genParticles;
	iEvent.getByLabel( genParticlesProducer_, genParticles );
	int iElectron=0; int iMuon=0; int iUnstableParticle=0;
	int iPartSel=0;  int iElectronSel=0; int iMuonSel=0; 
	int iJet=0, iMET=0, iJetSel=0, iMETSel=0;

	for(unsigned int j=0; j<genParticles->size(); ++j )
	{	
		const reco::GenParticle & p = (*genParticles)[ j ];
		//find the mother ID
		Int_t motherID = 0; Int_t grannyID = 0;
		if (p.numberOfMothers() > 0 )
		{
			//sanity check
			const Candidate * mom = p.mother();
			motherID = mom->pdgId();
			if (mom->numberOfMothers() > 0)
			{
				const Candidate * granny = mom->mother();
				grannyID = granny->pdgId();
				if ( motherID == p.pdgId() )
				{
					//check if the particle is "daugther of itself"
					motherID = granny->pdgId();
					if (granny->numberOfMothers() > 0) grannyID = (granny->mother())->pdgId();
				}
			}
		}
		
		if ( doElectronMC_ && abs(p.pdgId()) == 11 && p.status()==1 )
		{
			iElectron++;
			if ( abs(p.eta()>electronMC_etaMax_) || p.pt()<electronMC_ptMin_ ) continue;
			new( (*rootMCParticles)[iPartSel] ) cat::MCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(), p.vz(), p.pdgId(), p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, 0, 0, 0, 0, j );
			if(verbosity_>2) cout << "   ["<< setw(3) << iPartSel << "] MC Electron  " << (const cat::MCParticle&)(*rootMCParticles->At(iPartSel)) << endl;
			iPartSel++;
			iElectronSel++;
		}

		if ( doMuonMC_ && abs(p.pdgId()) == 13 && p.status()==1 )
		{
			iMuon++;
			if ( abs(p.eta()>muonMC_etaMax_) || p.pt()<muonMC_ptMin_ ) continue;
			new( (*rootMCParticles)[iPartSel] ) cat::MCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(), p.vz(), p.pdgId(), p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, 0, 0, 0, 0, j );
			if(verbosity_>2) cout << "   ["<< setw(3) << iPartSel << "] MC Muon  " << (const cat::MCParticle&)(*rootMCParticles->At(iPartSel)) << endl;
			iPartSel++;
			iMuonSel++;
		}

		// FIXME - GenJet collection instead
		if ( doJetMC_ && (abs(p.pdgId()) < 7 || abs(p.pdgId()) == 21 )&& p.status()==1 )
		{
			iJet++;
			if ( abs(p.eta()>jetMC_etaMax_) || p.pt()<jetMC_ptMin_ ) continue;
			new( (*rootMCParticles)[iPartSel] ) cat::MCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(), p.vz(), p.pdgId() , p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, 0, 0, 0, 0, j );
			if(verbosity_>2) cout << "   MC Jet  " << (const cat::Particle&)(*rootMCParticles->At(iPartSel)) << endl;
			iPartSel++;
			iJetSel++;
		}
    else if ( doJetMC_ && abs(p.pdgId()) == 5 && p.status() == 2 )
    {
      iJet++;
      if ( abs(p.eta()>jetMC_etaMax_) || p.pt()<jetMC_ptMin_ ) continue;
      new( (*rootMCParticles)[iPartSel] ) cat::MCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(), p.vz(), p.pdgId() , p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, 0, 0, 0, 0, j );
      if(verbosity_>2) cout << "   MC Jet  " << (const cat::Particle&)(*rootMCParticles->At(iPartSel)) << endl;
      iPartSel++;
      iJetSel++;
    }

		// FIXME - GenMET collection instead
		if ( doMETMC_ && (abs(p.pdgId()) == 12 || abs(p.pdgId()) == 14 ||  abs(p.pdgId()) == 16 || ( abs(p.pdgId()) > 1000000 && abs(p.pdgId()) < 3000000 ) )&& p.status()==1 )
		{
			iMET++;
			new( (*rootMCParticles)[iPartSel] ) cat::MCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(), p.vz(), p.pdgId() , p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, 0, 0, 0, 0, j );
			if(verbosity_>2) cout << "   MC MET  " << (const cat::Particle&)(*rootMCParticles->At(iPartSel)) << endl;
			iPartSel++;
			iMETSel++;
		}

		// add information on primary unstable particles: keep quarks, taus, Z, W, Higgs, susy and vlq particles, with status 3
		if ( doUnstablePartsMC_ && (abs(p.pdgId()) < 38 || (abs(p.pdgId()) > 1000000 && abs(p.pdgId()) < 3000000)  || (abs(p.pdgId()) > 4000000 && abs(p.pdgId()) < 6000000))	&& p.status()==3 )
		{
			iUnstableParticle++;	
			Int_t daug0Id = 0;
			Int_t daug1Id = 0;
			Int_t daug2Id = 0;
			Int_t daug3Id = 0;
			if (p.numberOfDaughters() > 0) daug0Id = p.daughter( 0 )->pdgId();
			if (p.numberOfDaughters() > 1) daug1Id = p.daughter( 1 )->pdgId();
			if (p.numberOfDaughters() > 2) daug2Id = p.daughter( 2 )->pdgId();
			if (p.numberOfDaughters() > 3) daug3Id = p.daughter( 3 )->pdgId();

			new( (*rootMCParticles)[iPartSel] ) cat::MCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(),p.vz(), p.pdgId(), p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, daug0Id, daug1Id, daug2Id, daug3Id, j );
			if(verbosity_>2) cout << "   ["<< setw(3) << iPartSel << "] unstable particle  " << (const cat::MCParticle&)(*rootMCParticles->At(iPartSel)) << endl;
			iPartSel++;		

		}	

	}
	
	if(verbosity_>1)
	{
		cout << endl;
		cout << "   Number of MC electrons = " << iElectron << ", preselected = " << iElectronSel << endl;
		cout << "   Number of MC muons = " << iMuon << ", preselected = " << iMuonSel << endl;
		cout << "   Number of primary unstable particles dumped in the ntuple = " << iUnstableParticle << endl;
	}	

}



