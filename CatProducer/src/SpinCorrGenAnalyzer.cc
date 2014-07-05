#include "../interface/SpinCorrGenAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

SpinCorrGenAnalyzer::SpinCorrGenAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0)
{
	genEventProducer_ = producersNames.getParameter<edm::InputTag>("genEventProducer");
}

SpinCorrGenAnalyzer::SpinCorrGenAnalyzer(const edm::ParameterSet& producersNames, int verbosity):verbosity_(verbosity)
{
	genEventProducer_ = producersNames.getParameter<edm::InputTag>("genEventProducer");
}

SpinCorrGenAnalyzer::SpinCorrGenAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	genEventProducer_ = producersNames.getParameter<edm::InputTag>("genEventProducer");
}

SpinCorrGenAnalyzer::~SpinCorrGenAnalyzer()
{
}


TLorentzVector  MyP4toTLV(reco::Particle::LorentzVector a){ return TLorentzVector(a.px(), a.py(), a.pz(), a.energy());}

void SpinCorrGenAnalyzer::Process(const edm::Event& iEvent, TClonesArray* rootSpinCorrGen){

	edm::Handle < TtGenEvent > genEvt;
        iEvent.getByLabel(genEventProducer_, genEvt);


	//	if(verbosity_>1) std::cout << "   SpinCorrGen  "  << "   Label: " << spinCorrGenProducer_.label() << "   Instance: " << spinCorrGenProducer_.instance() << std::endl;
	if(verbosity_>1) std::cout << "   SpinCorrGen  "  << "   Label: " << genEventProducer_.label() << "   Instance: " << genEventProducer_.instance() << std::endl;

	//my code
 
	double cosThetaTLHel= -9999;
	double cosThetaTBHel= -9999;
	double cosThetaTQHel= -9999;
	double cosPhi= -9999;
	double topsZMFMass =-9999;

	if(genEvt->isSemiLeptonic(WDecay::kMuon)){
	  //	  reco::Particle::LorentzVector topsZMF(genEvt.top()->p4()+genEvt.topBar()->p4());
	  reco::Particle::LorentzVector topsZMF(genEvt->leptonicDecayTop()->p4()+ genEvt->hadronicDecayTop()->p4());
	  
	  // boost particle 4-vectors to tt ZMF
	  if(genEvt->leptonicDecayTop() && genEvt->hadronicDecayTop() && genEvt->singleLepton() && genEvt->hadronicDecayB() && genEvt->hadronicDecayQuark()&& genEvt->hadronicDecayQuarkBar() ){
	    
	    TLorentzVector leptonicDecayTop;
	    TLorentzVector hadronicDecayTop;
	    TLorentzVector hadronicDecayB;
	    TLorentzVector hadronicDecayQuark;
	    TLorentzVector hadronicDecayQuarkBar;
	    TLorentzVector Muon;
	    
	    hadronicDecayTop      = MyP4toTLV(genEvt->hadronicDecayTop()->p4());		
	    leptonicDecayTop      = MyP4toTLV(genEvt->leptonicDecayTop()->p4());		
	    hadronicDecayB        = MyP4toTLV(genEvt->hadronicDecayB()->p4());
	    hadronicDecayQuark    = MyP4toTLV(genEvt->hadronicDecayQuark()->p4());
	    hadronicDecayQuarkBar = MyP4toTLV(genEvt->hadronicDecayQuarkBar()->p4());
	    Muon                  = MyP4toTLV(genEvt->singleLepton()->p4());
	    
	    
	    TLorentzVector topsZMF = leptonicDecayTop+ hadronicDecayTop;
	    
	    // boost particle 4-vectors to tt ZMF
	    TLorentzVector tLeptonicZMF  = leptonicDecayTop;
	    TLorentzVector tHadronicZMF  = hadronicDecayTop;
	    TLorentzVector lLeptonicZMF  = Muon;
	    TLorentzVector bHadronicZMF  = hadronicDecayB; 
	    TLorentzVector q1HadronicZMF = hadronicDecayQuark;
	    TLorentzVector q2HadronicZMF = hadronicDecayQuarkBar;
	    
	    tLeptonicZMF.Boost(- topsZMF.BoostVector()); 
	    tHadronicZMF.Boost(- topsZMF.BoostVector()); 
	    lLeptonicZMF.Boost(- topsZMF.BoostVector());  
	    bHadronicZMF.Boost(- topsZMF.BoostVector());  
	    q1HadronicZMF.Boost(- topsZMF.BoostVector()); 
	    q2HadronicZMF.Boost(- topsZMF.BoostVector()); 

	    //--------------------------------------------------------------------------------
	    // build spin basis unit vectors
	    TVector3 leptHelZMF(tLeptonicZMF.Vect().Unit());    
	    TVector3 hadrHelZMF(tHadronicZMF.Vect().Unit()); // = -leptHelZMF
	    TVector3 q1HadZMF(q1HadronicZMF.Vect().Unit());
	    TVector3 q2HadZMF(q2HadronicZMF.Vect().Unit());
	    TVector3 qHadZMF(q1HadronicZMF.E() < q2HadronicZMF.E() ? 
			     q1HadZMF :
			     q2HadZMF);// only lower energy quark used
      
	    // boost 4-vectors to t(bar) rest frames
	    TLorentzVector lLeptonicTRest  = lLeptonicZMF; //, tLeptonicZMF.BoostToCM()));
	    TLorentzVector bHadronicTRest  = bHadronicZMF; //, tHadronicZMF.BoostToCM()));
	    TLorentzVector q1HadronicTRest = q1HadronicZMF;//, tHadronicZMF.BoostToCM()));
	    TLorentzVector q2HadronicTRest = q2HadronicZMF;//, tHadronicZMF.BoostToCM()));
	    TLorentzVector qHadronicTRest;
	    
	    lLeptonicTRest.Boost(- tLeptonicZMF.BoostVector());  
	    bHadronicTRest.Boost(- tHadronicZMF.BoostVector());   
	    q1HadronicTRest.Boost(- tHadronicZMF.BoostVector());   
	    q2HadronicTRest.Boost(- tHadronicZMF.BoostVector());   
	    qHadronicTRest.Boost( -tHadronicZMF.BoostVector());   
	    
	    qHadronicTRest  = (q1HadronicTRest.E() < q2HadronicTRest.E() ? 
			       q1HadronicTRest :
			       q2HadronicTRest);// only lower energy quark used
	    
     
	    // extract particle directions in t(bar) rest frames
	    TVector3 lDirectionTRest =lLeptonicTRest.Vect().Unit();
	    TVector3 bDirectionTRest =bHadronicTRest.Vect().Unit();
	    TVector3 q1DirectionTRest=q1HadronicTRest.Vect().Unit();
	    TVector3 q2DirectionTRest=q2HadronicTRest.Vect().Unit();
	    TVector3 qDirectionTRest =qHadronicTRest.Vect().Unit();
	    
	    cosThetaTLHel = leptHelZMF.Dot(lDirectionTRest);
	    cosThetaTBHel = hadrHelZMF.Dot(bDirectionTRest);
	    cosThetaTQHel = hadrHelZMF.Dot(qDirectionTRest);
	    cosPhi = lDirectionTRest.Dot(qDirectionTRest);
	    topsZMFMass = topsZMF.M();
	    
	    
	
	    //set the variable in the cat
	    cat::SpinCorrGen localSSG;
	    localSSG.setcosThetaTLHel(cosThetaTLHel); 
	    localSSG.setcosThetaTBHel(cosThetaTBHel);
	    localSSG.setcosThetaTQHel(cosThetaTQHel);
	    localSSG.setcosPhi(cosPhi);
	    localSSG.settopsZMFMass(topsZMFMass);
	    
	    new( (*rootSpinCorrGen)[0] ) cat::SpinCorrGen(localSSG);
	    
	  }// close if 
	}//close if semilep


}
