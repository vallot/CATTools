#ifndef CATTools_DataFormats_NPGenEvent_h
#define CATTools_DataFormats_NPGenEvent_h
#include <iostream>
#include <iomanip>
#include <vector>

#include "TObject.h"
#include "Rtypes.h"

#include "../interface/MCParticle.h"
#include "../interface/GenTop.h"

namespace cat
{
	class NPGenEvent: public TObject{

	public:
		NPGenEvent(){isNewPhysics_ = false;};
		NPGenEvent(const Bool_t isNewPhysics, const std::vector<cat::GenTop> tops, const std::vector<cat::MCParticle> leptons, const std::vector<cat::MCParticle> quarks,
                    const std::vector<cat::MCParticle> bquarks, const std::vector<cat::MCParticle> invisibleParticles, const std::vector<cat::MCParticle> neutrinos,
                    const std::vector<cat::MCParticle> gluinos, const std::vector<cat::MCParticle> stops)
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
		NPGenEvent(const NPGenEvent& gevt)
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

		virtual ~NPGenEvent(){;};

		virtual TString typeName() const { return "NPGenEvent"; }

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

		std::vector<GenTop> tops() const {return tops_;};
		std::vector<MCParticle> leptons() const {return leptons_;};
		std::vector<MCParticle> quarks() const {return quarks_;};
		std::vector<MCParticle> bquarks() const {return bquarks_;};
		std::vector<MCParticle> invisibleParticles() const {return invisibleParticles_;};
		std::vector<MCParticle> neutrinos() const {return neutrinos_;}; 
		std::vector<MCParticle> gluinos() const {return gluinos_;}; 
		std::vector<MCParticle> stops() const {return stops_;}; 

		friend std::ostream& operator<< (std::ostream& stream, const NPGenEvent& gevent)
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
    std::vector<GenTop> tops_;
		std::vector<MCParticle> leptons_;
		std::vector<MCParticle> quarks_;
		std::vector<MCParticle> bquarks_;
		std::vector<MCParticle> invisibleParticles_;
		std::vector<MCParticle> neutrinos_;
		std::vector<MCParticle> gluinos_;
		std::vector<MCParticle> stops_;
  
		ClassDef (NPGenEvent,3); 
	};
}

#endif
