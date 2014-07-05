#ifndef Event_h
#define Event_h

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "Rtypes.h"
#include "TObject.h"
#include "TVector3.h"

using namespace std;

namespace cat
{

	struct triggeredObject{
		triggeredObject() {};
                virtual ~triggeredObject(){};
		Int_t id;
		Float_t pt;
		Float_t eta;
		Float_t phi;

		ClassDef(triggeredObject,2);

	};

	
	class Event : public TObject
	{
	
	public:
	
		Event() : 
			nb_(0)
			,eventId_(-1)
			,runId_(-1)
			,lumiBlockId_(-1)
			,flavHistPath_(-1)
			,kt6PFJets_rho_(-1)
      			,kt6PFJetsForIsolation_rho_(-1)
			,nTracks_(-1)
			,nHighPurityTracks_(-1)
			,nTruePU_(-1.)
			,passGlobalHLT_(false)
			,trigHLT_(0)
			,idParton1_(-1)
			,xParton1_(-1.)
			,idParton2_(-1)
			,xParton2_(-1.)
			,factorizationScale_(-1.)
			{;}

		~Event() {;}
	
		// Event number
		Int_t nb() const { return nb_; }
		Int_t eventId() const { return eventId_; }
		Int_t runId() const { return runId_; }
		Int_t lumiBlockId() const {return lumiBlockId_; }
		Int_t flavorHistoryPath() const { return flavHistPath_; }
		double kt6PFJets_rho() const { return kt6PFJets_rho_; }
    		double kt6PFJetsForIsolation_rho() const { return kt6PFJetsForIsolation_rho_; }
		// generalTracks infos (for cleaning of the scraping events)
		Int_t nTracks() const { return nTracks_; }
		Int_t nHighPurityTracks() const { return nHighPurityTracks_; }
		// PileUp info
		Int_t nPu(Int_t nBX) { return ( nPu_.find(nBX) != nPu_.end() ? nPu_[nBX] : -9999 ); }
		Float_t nTruePU() { return nTruePU_; }
		// Trigger decision
		Bool_t passGlobalHLT() const { return passGlobalHLT_; }
		UInt_t nHLTPaths() const { return trigHLT_.size(); }
		std::vector<Bool_t> trigHLT() const { return trigHLT_; }
		Bool_t trigHLT(UInt_t i) const
		{
			if (trigHLT_.size()>i)
			{
				return trigHLT_.at(i); 
			}
			else
			{
				cout << "HLT path " << i << " not found" << endl;
				return false;
			}
		}
		
		// PDF infos
		// flavour first incoming parton
		Int_t idParton1() const { return idParton1_; }
		// energy fraction carried by first incoming parton
		Float_t xParton1() const { return xParton1_; }
		// flavour first incoming parton
		Int_t idParton2() const { return idParton2_; }
		// energy fraction carried by first incoming parton
		Float_t xParton2() const { return xParton2_; }
		// Q-scale used in evaluation of PDF's (in GeV). 
		Float_t factorizationScale() const { return factorizationScale_; }

		std::map<std::string, std::vector<cat::triggeredObject> >  getTriggerFilters() const { return triggerFilters_;}


		void setNb(Int_t nb) { nb_ = nb; }
		void setEventId(Int_t eventId) { eventId_ = eventId; }
		void setRunId(Int_t runId) { runId_ = runId; }
		void setLumiBlockId(Int_t lumiBlockId) { lumiBlockId_ = lumiBlockId; }
		void setflavorHistoryPath(Int_t flavHistPath) { flavHistPath_ = flavHistPath; }
		void setKt6PFJets_rho(double kt6PFJets_rho) { kt6PFJets_rho_ = kt6PFJets_rho; }
		void setKt6PFJetsForIsolation_rho(double kt6PFJetsForIsolation_rho) { kt6PFJetsForIsolation_rho_ = kt6PFJetsForIsolation_rho; }
		void setNTracks(Int_t nTracks) { nTracks_ = nTracks; }
		void setNPu(Int_t nBX, Int_t nPu) { nPu_[nBX] = nPu; }
		void setNTruePU(Float_t nTruePU) { nTruePU_ = nTruePU; }
		void setNHighPurityTracks(Int_t nHighPurityTracks) { nHighPurityTracks_ = nHighPurityTracks; }

		void setGlobalHLT(Bool_t passGlobalHLT) { passGlobalHLT_ = passGlobalHLT; }
		void setTrigHLT(std::vector<Bool_t> trigHLT)
		{
			trigHLT_.resize(trigHLT.size());
			for (unsigned int i=0; i!=trigHLT.size(); ++i) trigHLT_[i]=trigHLT[i];
		}

		void setIdParton1(Int_t idParton1) { idParton1_=idParton1; }
		void setXParton1(Float_t xParton1) { xParton1_=xParton1; }
		void setIdParton2(Int_t idParton2) { idParton2_=idParton2; }
		void setXParton2(Float_t xParton2) { xParton2_=xParton2; }
		void setFactorizationScale(Float_t factorizationScale) { factorizationScale_=factorizationScale; }

		void AddTriggerObject(string path, Int_t id, Float_t pt, Float_t eta, Float_t phi){
			cat::triggeredObject object;
			object.id = id;
			object.pt = pt;
			object.eta = eta;
			object.phi = phi;
			triggerFilters_[path].push_back(object);
		}
		
		/*
		// FIXME
		friend std::ostream& operator<< (std::ostream& stream, const Event& event) {
			stream << "Event #"<< event.nb() <<"  HLT="<< event.trigHLT();
			return stream;
		};
		*/
	
	private:
	
		Int_t nb_;
		Int_t eventId_;
		Int_t runId_;
		Int_t lumiBlockId_;
		
		Int_t flavHistPath_; // flavor history path, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFlavorHistory
		double kt6PFJets_rho_;
		double kt6PFJetsForIsolation_rho_;
		
		// generalTracks infos (for cleaning of the scraping events)
		Int_t nTracks_;
		Int_t nHighPurityTracks_;
		
		// PileUp info
		std::map<Int_t,Int_t> nPu_;
		Float_t nTruePU_;
		// Trigger Infos
		Bool_t passGlobalHLT_;
		std::vector<Bool_t> trigHLT_;

		// PDF infos
		Int_t idParton1_;
		Float_t xParton1_;
		Int_t idParton2_;
		Float_t xParton2_;
		Float_t factorizationScale_;

		std::map<std::string, std::vector<cat::triggeredObject> > triggerFilters_;

		ClassDef (Event,4);
	};
}

#endif
