#ifndef CatNPGenEvent_h
#define CatNPGenEvent_h
#include <iostream>
#include <iomanip>
#include <vector>

#include "TObject.h"
#include "Rtypes.h"

#include "../interface/CatMCParticle.h"
#include "../interface/CatGenTop.h"

namespace cat
{
	class CatNPGenEvent: public TObject{

	public:
		CatNPGenEvent(){isNewPhysics_ = false;};
		CatNPGenEvent(const Bool_t isNewPhysics, const std::vector<cat::CatGenTop> tops, const std::vector<cat::CatMCParticle> leptons, const std::vector<cat::CatMCParticle> quarks,
                    const std::vector<cat::CatMCParticle> bquarks, const std::vector<cat::CatMCParticle> invisibleParticles, const std::vector<cat::CatMCParticle> neutrinos,
                    const std::vector<cat::CatMCParticle> gluinos, const std::vector<cat::CatMCParticle> stops)
		{
			isNewPhysics_ = isNewPhysics;
			tops_ = tops;
			leptons_ = leptons;
			quarks_ = quarks;
			bquarks_ = bquarks;
			invisibleParticles_ = invisibleParticles;
			neutrinos_ = neutrinos;
			gluinos_ = gluinos;
			stops_ = stops;
		};
		CatNPGenEvent(const CatNPGenEvent& gevt)
		{
			isNewPhysics_ = gevt.isNewPhysics_;
			tops_ = gevt.tops_;
			leptons_ = gevt.leptons_;
			quarks_ = gevt.quarks_;
			bquarks_ = gevt.bquarks_;
			invisibleParticles_ = gevt.invisibleParticles_;
			neutrinos_ = gevt.neutrinos_;
			gluinos_ = gevt.gluinos_;
			stops_ = gevt.stops_;
		};

		virtual ~CatNPGenEvent(){;};

		virtual TString typeName() const { return "CatNPGenEvent"; }

		Bool_t isNewPhysics() const {return isNewPhysics_;};
		Bool_t isThereTop() const {return tops_.size()>0? true:false;};

		Int_t numberOfTops() const {return tops_.size();};
		Int_t numberOfLeptons() const {return leptons_.size();};
		Int_t numberOfQuarks() const {return quarks_.size();};
		Int_t numberOfBQuarks() const {return bquarks_.size();};
		Int_t numberOfInvisibleParticles() const {return invisibleParticles_.size();};
		Int_t numberOfNeutrinos() const {return neutrinos_.size();};
		Int_t numberOfGluinos() const{ return gluinos_.size();};
		Int_t numberOfStops() const{ return stops_.size();};

		std::vector<CatGenTop> tops() const {return tops_;};
		std::vector<CatMCParticle> leptons() const {return leptons_;};
		std::vector<CatMCParticle> quarks() const {return quarks_;};
		std::vector<CatMCParticle> bquarks() const {return bquarks_;};
		std::vector<CatMCParticle> invisibleParticles() const {return invisibleParticles_;};
		std::vector<CatMCParticle> neutrinos() const {return neutrinos_;}; 
		std::vector<CatMCParticle> gluinos() const {return gluinos_;}; 
		std::vector<CatMCParticle> stops() const {return stops_;}; 

		friend std::ostream& operator<< (std::ostream& stream, const CatNPGenEvent& gevent)
		{
			stream <<"Event ";
			if (gevent.isNewPhysics()) stream << " is NewPhysics ";
			stream <<". ";
			stream << gevent.numberOfTops() <<" tops - "<< gevent.numberOfLeptons()<< " leptons - "<< gevent.numberOfQuarks() <<" quarks - ( "<< gevent.numberOfBQuarks()<<"b ) ";
			stream << gevent.numberOfInvisibleParticles() <<" invisible particles ( "<< gevent.numberOfNeutrinos() <<" neutrinos )";
			return stream;
		}

	private:
  
		Bool_t isNewPhysics_;
    std::vector<CatGenTop> tops_;
		std::vector<CatMCParticle> leptons_;
		std::vector<CatMCParticle> quarks_;
		std::vector<CatMCParticle> bquarks_;
		std::vector<CatMCParticle> invisibleParticles_;
		std::vector<CatMCParticle> neutrinos_;
		std::vector<CatMCParticle> gluinos_;
		std::vector<CatMCParticle> stops_;
  
		ClassDef (CatNPGenEvent,3); 
	};
}

#endif
