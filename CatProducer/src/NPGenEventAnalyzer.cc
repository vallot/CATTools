#include "../interface/NPGenEventAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

NPGenEventAnalyzer::NPGenEventAnalyzer (const edm::ParameterSet & producersNames):
verbosity_ (0)
{
  genParticlesProducer_ = producersNames.getParameter < edm::InputTag > ("genParticlesProducer");
}

NPGenEventAnalyzer::NPGenEventAnalyzer (const edm::ParameterSet & producersNames, int verbosity):
verbosity_ (verbosity)
{
  genParticlesProducer_ = producersNames.getParameter < edm::InputTag > ("genParticlesProducer");
}

NPGenEventAnalyzer::NPGenEventAnalyzer (const edm::ParameterSet & producersNames, const edm::ParameterSet & myConfig, int verbosity):
verbosity_ (verbosity)
{
  genParticlesProducer_ = producersNames.getParameter < edm::InputTag > ("genParticlesProducer");
}

NPGenEventAnalyzer::~NPGenEventAnalyzer ()
{
}

TLorentzVector NPGenEventAnalyzer::P4toTLV (reco::Particle::LorentzVector a)
{
  return TLorentzVector (a.px (), a.py (), a.pz (), a.energy ());
}

CatMCParticle NPGenEventAnalyzer::ConvertMCPart(reco::GenParticleCollection::const_iterator t){
  TLorentzVector p4 = P4toTLV(t->p4());
  TVector3 vertex = TVector3(t->vertex().x(),t->vertex().y(),t->vertex().z());
  int motherType = -9999;
  int grannyType = -9999;
  if(t->mother()) motherType = t->mother()->pdgId();
  if(t->mother() && t->mother()->mother()) grannyType = t->mother()->mother()->pdgId();
  int dauOneId = -9999;
  int dauTwoId = -9999;
  int dauThreeId = -9999;
  int dauFourId = -9999;
  int nDau = 0 ;
  for (reco::GenParticle::const_iterator td = t->begin (); td != t->end (); ++td){
    nDau++;
    if(nDau==1) dauOneId = td->pdgId();
    if(nDau==2) dauTwoId = td->pdgId();
    if(nDau==3) dauThreeId = td->pdgId();
    if(nDau==4) dauFourId = td->pdgId();
  }
  CatMCParticle part(p4, vertex, t->pdgId(), t->charge(), t->status(), nDau, motherType, grannyType, dauOneId, dauTwoId, dauThreeId, dauFourId, -9999);
  return part;
}

CatMCParticle NPGenEventAnalyzer::ConvertMCPart(reco::GenParticle::const_iterator t){
  TLorentzVector p4 = P4toTLV(t->p4());
  TVector3 vertex = TVector3(t->vertex().x(),t->vertex().y(),t->vertex().z());
  int motherType = -9999;
  int grannyType = -9999;
  if(t->mother()) motherType = t->mother()->pdgId();
  if(t->mother() && t->mother()->mother()) grannyType = t->mother()->mother()->pdgId();
  int dauOneId = -9999;
  int dauTwoId = -9999;
  int dauThreeId = -9999;
  int dauFourId = -9999;
  int nDau = 0 ;
  for (reco::GenParticle::const_iterator td = t->begin (); td != t->end (); ++td){
    nDau++;
    if(nDau==1) dauOneId = td->pdgId();
    if(nDau==2) dauTwoId = td->pdgId();
    if(nDau==3) dauThreeId = td->pdgId();
    if(nDau==4) dauFourId = td->pdgId();
  }
  CatMCParticle part(p4, vertex, t->pdgId(), t->charge(), t->status(), nDau, motherType, grannyType, dauOneId, dauTwoId, dauThreeId, dauFourId, -9999);
  return part;
}
void
NPGenEventAnalyzer::Process (const edm::Event & iEvent, TClonesArray * rootGenEvent)
{

  //cout<<"Handle"<<endl;
  edm::Handle < reco::GenParticleCollection > hsrc;
  //cout<<"Handle"<<endl;
  iEvent.getByLabel (genParticlesProducer_, hsrc);
  //cout<<"Handle"<<endl;
  reco::GenParticleCollection src = * hsrc;
  //cout<<"Handle"<<endl;

  if (verbosity_ > 1)
    std::cout << "   NPGenEventAnalyzer  " << "   Label: " << genParticlesProducer_.label () << "   Instance: " << genParticlesProducer_.instance () << std::endl;

  bool isNewPhysics_ = false;
  vector < CatGenTop > tops_;
  vector < CatMCParticle > stops_;
  vector < CatMCParticle > leptons_;
  vector < CatMCParticle > quarks_;
  vector < CatMCParticle > bquarks_;
  vector < CatMCParticle > gluinos_;
  vector < CatMCParticle > neutrinos_;
  vector < CatMCParticle > invisibleParticles_;

  for (reco::GenParticleCollection::const_iterator t = src.begin (); t != src.end (); ++t)
    {
      if(verbosity_>4) cout<<t->p4()<<" "<<t->pdgId()<<" "<<t->status()<<endl;
      if (t->status () == 3)
	{
	  if (t->pdgId () > 100000)
	    {
	      isNewPhysics_ = true;
	    }
	  if ((abs (t->pdgId ()) == 11 || abs (t->pdgId ()) == 13 || abs (t->pdgId ()) == 15) &&
	      ((abs (t->mother ()->pdgId ()) > 22 && abs (t->mother ()->pdgId ()) < 40) || (abs (t->mother ()->pdgId ()) > 999999 && abs (t->mother ()->pdgId ()) < 2999999)))
	    leptons_.push_back (ConvertMCPart(t));
	  if ((abs (t->pdgId ()) == 12 || abs (t->pdgId ()) == 14 || abs (t->pdgId ()) == 16) &&
	      ((abs (t->mother ()->pdgId ()) > 22 && abs (t->mother ()->pdgId ()) < 40) || (abs (t->mother ()->pdgId ()) > 999999 && abs (t->mother ()->pdgId ()) < 2999999)))
	    neutrinos_.push_back (ConvertMCPart(t));
	  if (abs (abs (t->pdgId ())) == 1000022 || abs (abs (t->pdgId ())) == 1000023 || abs (abs (t->pdgId ())) == 1000025 || abs (abs (t->pdgId ())) == 1000035)
	    {
	      int nofd = 0;
	      for (reco::GenParticle::const_iterator td = t->begin (); td != t->end (); ++td)
		nofd++;
	      if (nofd < 2)
		invisibleParticles_.push_back (ConvertMCPart(t));
	    }
	  if (abs (t->pdgId ()) < 6)
	    quarks_.push_back (ConvertMCPart(t));
	  if (abs (t->pdgId ()) == 5)
	    bquarks_.push_back (ConvertMCPart(t));
	  if (abs (t->pdgId ()) == 1000021)
	    gluinos_.push_back (ConvertMCPart(t));
	  if (abs (t->pdgId ()) == 1000006 || abs (t->pdgId ()) == 2000006)
	    stops_.push_back (ConvertMCPart(t));
	  if (abs (t->pdgId ()) == 6)
	    {
	      bool isLeptonic_ = false;
	      CatMCParticle top_;
	      CatMCParticle W_;
	      CatMCParticle bquark_;
	      CatMCParticle quark_;
	      CatMCParticle quarkBar_;
	      CatMCParticle lepton_;
	      CatMCParticle neutrino_;
	      
	      top_ = ConvertMCPart(t);
	      reco::GenParticle::const_iterator td = t->begin ();
	      for (; td != t->end (); ++td)
		{
		  if (abs (td->pdgId ()) == 24)
		    {
		      W_  = ConvertMCPart(td);
		      GenParticle::const_iterator Wd = td->begin ();
		      for (; Wd != td->end (); ++Wd)
			{
			  if (Wd->pdgId () > 0 && Wd->pdgId () < 6)
			    {
			      quark_  = ConvertMCPart(Wd);
			    }
			  if (Wd->pdgId () < 0 && Wd->pdgId () > -6)
			    {
			      quarkBar_  = ConvertMCPart(Wd);
			    }
			  if (abs (Wd->pdgId ()) == 11 || abs (Wd->pdgId ()) == 13 || abs (Wd->pdgId ()) == 15)
			    {
			      lepton_  = ConvertMCPart(Wd);
			      isLeptonic_ = true;
			    }
			  if (abs (Wd->pdgId ()) == 12 || abs (Wd->pdgId ()) == 14 || abs (Wd->pdgId ()) == 16)
			    {
			      neutrino_  = ConvertMCPart(Wd);
			    }
			}
		    }
		  else
		    bquark_  = ConvertMCPart(td);
		}
	      //cout<<"TopGenPart"<<endl;
	      //
	      char cprod[1000];
	      char cprodTemp[1000];
	      sprintf(cprod,"%d->",t->motherRef()->pdgId());
	      GenParticle::const_iterator tm = t->motherRef()->begin();
	      for( ; tm!=t->motherRef()->end(); ++tm){
	         sprintf(cprodTemp,"%s %d",cprod,tm->pdgId());
		 sprintf(cprod,"%s",cprodTemp);
	      }
	      string production_ = string(cprod);	  
	      if(isLeptonic_){
	        CatGenTop top (isLeptonic_, top_, W_, bquark_, lepton_, neutrino_, production_);
	        tops_.push_back (top);
	      }
	      else{
	        CatGenTop top (isLeptonic_, top_, W_, bquark_, quark_, quarkBar_, production_);
	        tops_.push_back (top);
	     }
	      //cout<<"End TopGenPart"<<endl;
	    }
	}
    }
    //cout<<"AT THE END"<<endl;
  CatNPGenEvent CatgenEvt (isNewPhysics_, tops_, leptons_, quarks_, bquarks_, invisibleParticles_, neutrinos_, gluinos_, stops_);
    //cout<<"AT THE END"<<endl;


  new ((*rootGenEvent)[0]) CatNPGenEvent (CatgenEvt);
    //cout<<"AT THE END"<<endl;
}
