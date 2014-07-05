#include "../interface/GenEventAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

GenEventAnalyzer::GenEventAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	genEventProducer_ = producersNames.getParameter<edm::InputTag>("genEventProducer");
}

GenEventAnalyzer::GenEventAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	genEventProducer_ = producersNames.getParameter<edm::InputTag>("genEventProducer");
}

GenEventAnalyzer::GenEventAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	genEventProducer_ = producersNames.getParameter<edm::InputTag>("genEventProducer");
}

GenEventAnalyzer::~GenEventAnalyzer()
{
}


TLorentzVector  P4toTLV(reco::Particle::LorentzVector a){ return TLorentzVector(a.px(), a.py(), a.pz(), a.energy());}

void GenEventAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootGenEvent){

	edm::Handle < TtGenEvent > genEvent;
        iEvent.getByLabel(genEventProducer_, genEvent);


	if(verbosity_>1) std::cout << "   GenEvent  "  << "   Label: " << genEventProducer_.label() << "   Instance: " << genEventProducer_.instance() << std::endl;

        cat::CatGenEvent genEvt;
	genEvt.SetBoolean(genEvent->isTtBar(), genEvent->isFullHadronic(), genEvent->isSemiLeptonic(), genEvent->isFullLeptonic());
	genEvt.SetSemiLeptonicChannel((cat::CatGenEvent::LeptonType) genEvent->semiLeptonicChannel());
	TLorentzVector neutrino;
	TLorentzVector lepton;
	TLorentzVector leptonicDecayW, leptonicDecayB, leptonicDecayTop;
	TLorentzVector hadronicDecayW, hadronicDecayB, hadronicDecayTop, hadronicDecayQuark, hadronicDecayQuarkBar;
	if( genEvent->isSemiLeptonic()){
	  if(genEvent->singleNeutrino()) neutrino = P4toTLV(genEvent->singleNeutrino()->p4());
  	  if(genEvent->singleLepton()) lepton = P4toTLV(genEvent->singleLepton()->p4());
	  if(genEvent->leptonicDecayW() && genEvent->wMinus() && genEvent->wPlus()) leptonicDecayW = P4toTLV(genEvent->leptonicDecayW()->p4());
	  if(genEvent->leptonicDecayB()) leptonicDecayB = P4toTLV(genEvent->leptonicDecayB()->p4());
	  if(genEvent->leptonicDecayTop()) leptonicDecayTop = P4toTLV(genEvent->leptonicDecayTop()->p4());
	  //cerr << "Adress of Wminus " << genEvent->wMinus() << " " << genEvent->wMinus()->p4() << endl;
	  //cerr << "Adress of Wplus " << genEvent->wPlus() << " " << genEvent->wPlus()->p4() << endl;

	  if(genEvent->hadronicDecayW() && genEvent->wMinus() && genEvent->wPlus()) hadronicDecayW = P4toTLV(genEvent->hadronicDecayW()->p4()); 
	  if(genEvent->singleLepton() && genEvent->singleLepton()->pdgId()>0 && genEvent->wMinus()) hadronicDecayW = P4toTLV(genEvent->wMinus()->p4());
 	  if(genEvent->singleLepton() && genEvent->singleLepton()->pdgId()<0 && genEvent->wPlus()) hadronicDecayW = P4toTLV(genEvent->wPlus()->p4());

	  if(genEvent->hadronicDecayB()) hadronicDecayB = P4toTLV(genEvent->hadronicDecayB()->p4());
	  if(genEvent->hadronicDecayTop()) hadronicDecayTop = P4toTLV(genEvent->hadronicDecayTop()->p4());
	  if(genEvent->hadronicDecayQuark()) hadronicDecayQuark = P4toTLV(genEvent->hadronicDecayQuark()->p4());
	  if(genEvent->hadronicDecayQuarkBar()) hadronicDecayQuarkBar = P4toTLV(genEvent->hadronicDecayQuarkBar()->p4());
	}
	genEvt.SetTLorentzVector(lepton, neutrino, leptonicDecayW, leptonicDecayB, leptonicDecayTop, hadronicDecayW, hadronicDecayB, hadronicDecayTop, hadronicDecayQuark, hadronicDecayQuarkBar);
	std::vector<TLorentzVector> ISR, leptonicDecayTopRadiation, hadronicDecayTopRadiation;
	for(unsigned int i=0;i<genEvent->topSisters().size();i++) 
	 if(genEvent->topSisters()[i]) 
	   ISR.push_back(P4toTLV(genEvent->topSisters()[i]->p4())); 
	if( genEvent->isSemiLeptonic()){
	 for(unsigned int i=0;i<genEvent->radiatedGluons(genEvent->leptonicDecayTop()->pdgId()).size();i++) 
           if(genEvent->radiatedGluons(genEvent->leptonicDecayTop()->pdgId())[i]) 
	    leptonicDecayTopRadiation.push_back(P4toTLV(genEvent->radiatedGluons(genEvent->leptonicDecayTop()->pdgId())[i]->p4())); 
	 for(unsigned int i=0;i<genEvent->radiatedGluons(genEvent->hadronicDecayTop()->pdgId()).size();i++) 
	  if(genEvent->radiatedGluons(genEvent->hadronicDecayTop()->pdgId())[i])
	    hadronicDecayTopRadiation.push_back(P4toTLV(genEvent->radiatedGluons(genEvent->hadronicDecayTop()->pdgId())[i]->p4())); 
	}
	genEvt.SetRadiation(leptonicDecayTopRadiation, hadronicDecayTopRadiation, ISR);
	new( (*rootGenEvent)[0] ) cat::CatGenEvent(genEvt);
}
