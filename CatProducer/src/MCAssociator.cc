#include "../interface/MCAssociator.h"

using namespace cat;

MCAssociator::MCAssociator(): verbosity_(0), nMC_(0), mcParticles_(0), genParticles_(), mcParticlesMap_()
{
}


MCAssociator::MCAssociator(const edm::ParameterSet& producersNames, int verbosity) : verbosity_(verbosity), nMC_(0), mcParticles_(0), genParticles_(), mcParticlesMap_()
{
	genParticlesProducer_ = producersNames.getParameter<edm::InputTag>("genParticlesProducer");
}



void MCAssociator::init(const edm::Event& iEvent, TClonesArray* mcParticles)
{
	// FIXME - Protect against no genParticles collection in PoolSource
	iEvent.getByLabel( genParticlesProducer_, genParticles_ );

	// fill map<igen,imc> where igen=index in genParticle collection and imc=index in mcParticles TClonesArray
	if(verbosity_>1) cout << endl << "Matching recoParticles to mcParticles... " << endl;
	int igen;
	mcParticles_ = mcParticles;
	nMC_ = mcParticles_->GetEntriesFast();
	for (int imc=0; imc<nMC_; imc++)
	{
		// TODO - remove indexInList in MCParticle
		//igen = ( (cat::Particle*)mcParticles_->At(imc))->genParticleIndex_();
		igen = ( (cat::MCParticle*)mcParticles_->At(imc))->genParticleIndex();
		mcParticlesMap_[igen]=imc;
	}
}


void MCAssociator::process(TClonesArray* recoParticles)
{
	int igen=-1;
	std::map<int,int>::iterator it;

	for (int ipart=0; ipart<recoParticles->GetEntriesFast(); ipart++)
	{
		cat::Particle* recoParticle = (cat::Particle*)recoParticles->At(ipart);
		igen = recoParticle->genParticleIndex();
		if(igen<=0)
		{
			if(verbosity_>2) cout <<"   "<< recoParticle->typeName() << "[" << ipart << "] not matched (at CMSSW level)" << endl;
			continue;
		}
		
		it=mcParticlesMap_.find(igen);
		if(it==mcParticlesMap_.end())
		{
			if(verbosity_>2) cout <<"   "<< recoParticle->typeName() << "[" << ipart << "] not matched (at TotoAna level)... add new cat::MCParticle..." << endl;
			// if igen not found in mcParticles[], add genParticle in mcParticles[]
			const reco::GenParticle & p = (*genParticles_)[igen];
			//find the mother ID
			Int_t motherID = 0; Int_t grannyID = 0;
			if (p.numberOfMothers() > 0 )
			{
				//sanity check
				const reco::Candidate * mom = p.mother();
				motherID = mom->pdgId();
				if (mom->numberOfMothers() > 0)
				{
					const reco::Candidate * granny = mom->mother();
					grannyID = granny->pdgId();
					if ( motherID == p.pdgId() )
					{
						//check if the particle is "daugther of itself"
						motherID = granny->pdgId();
						if (granny->numberOfMothers() > 0) grannyID = (granny->mother())->pdgId();
					}
				}
			}

			cat::MCParticle localMCParticle( p.px(), p.py(), p.pz(), p.energy(), p.vx(), p.vy(), p.vz(), p.pdgId(), p.charge(), p.status(), p.numberOfDaughters(), motherID, grannyID, 0, 0, 0, 0, igen );
			new( (*mcParticles_)[nMC_] ) cat::MCParticle(localMCParticle);
			if(verbosity_>2) cout <<"      ===> now matched to mcParticle["<< nMC_<<"] "<< localMCParticle << endl;
			nMC_++;
		}
		else
		{
			if(verbosity_>2) cout <<"   "<< ( (cat::Particle*)recoParticles->At(ipart))->typeName() << "[" << ipart << "] matched to mcParticles[" << it->second << "]" << endl;
		}
	}
}


void MCAssociator::printParticleAssociation(TClonesArray* recoParticles)
{
	for (Int_t imc=0; imc<nMC_; imc++) // reindexing to take into account matches that were not inserted with MCAssociator::Init
	{
		Int_t igen = ( (cat::MCParticle*)mcParticles_->At(imc))->genParticleIndex();
		mcParticlesMap_[igen]=imc;
	}
	std::map<int,int>::iterator it;
	for (int ipart=0; ipart<recoParticles->GetEntriesFast(); ipart++)
	{
		if (ipart==0) cout << endl;
		cat::Particle* localParticle = (cat::Particle*)recoParticles->At(ipart);
    Int_t indexRECO = localParticle->genParticleIndex();
		it=mcParticlesMap_.find(indexRECO);
		if(it!=mcParticlesMap_.end())
		{
			cout <<"   "<< localParticle->typeName() <<"["<< ipart << "] " << *localParticle<< endl;
			cout <<"       ===> matched to mcParticles: " << *((cat::MCParticle*)mcParticles_->At(it->second)) << endl;

		} else {
			cout <<"   "<< localParticle->typeName() <<"["<< ipart << "] " << *localParticle<< endl;
			cout <<"       ===> not matched" << endl;
		}
	}
}
