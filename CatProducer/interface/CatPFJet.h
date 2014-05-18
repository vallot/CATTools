#ifndef CatPFJet_h
#define CatPFJet_h

#include "../interface/CatParticle.h"
#include "../interface/CatJet.h"

#include "Rtypes.h"
#include "TObject.h"

#include <iostream>
#include <iomanip>

using namespace std;

namespace cat
{
	class CatPFJet : public CatJet
	{

	public:
		CatPFJet() :
			CatJet()
			,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		CatPFJet(const CatPFJet& jet) :
			CatJet(jet)
			,chargedHadronEnergyFraction_(jet.chargedHadronEnergyFraction_)
			,neutralHadronEnergyFraction_(jet.neutralHadronEnergyFraction_)
			,chargedEmEnergyFraction_(jet.chargedEmEnergyFraction_)
			,chargedMuEnergyFraction_(jet.chargedMuEnergyFraction_)
			,neutralEmEnergyFraction_(jet.neutralEmEnergyFraction_)			
			,HFHadronEnergyFraction_(jet.HFHadronEnergyFraction_)
			,HFEMEnergyFraction_(jet.HFEMEnergyFraction_)
			,chargedMultiplicity_(jet.chargedMultiplicity_)
			,neutralMultiplicity_(jet.neutralMultiplicity_)
			,muonMultiplicity_(jet.muonMultiplicity_)
			,HFHadronMultiplicity_(jet.HFHadronMultiplicity_)
			,HFEMMultiplicity_(jet.HFEMMultiplicity_)			
			{;}

		CatPFJet(const CatJet& jet) :
			CatJet(jet)		        
			,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		CatPFJet(Double_t px, Double_t py, Double_t pz, Double_t e) :
			CatJet(px,py,px,e)
		  ,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		CatPFJet(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			CatJet(px,py,pz,e,vtx_x,vtx_y,vtx_z)
		  ,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		CatPFJet(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
			CatJet(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
		  ,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		CatPFJet(const TLorentzVector &momentum) :
			CatJet(momentum)
			,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		CatPFJet(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
			CatJet(momentum, vertex, type, charge)
			,chargedHadronEnergyFraction_(-9999.)
			,neutralHadronEnergyFraction_(-9999.)
			,chargedEmEnergyFraction_(-9999.)
			,chargedMuEnergyFraction_(-9999.)
			,neutralEmEnergyFraction_(-9999.)			
			,HFHadronEnergyFraction_(-9999.)
			,HFEMEnergyFraction_(-9999.)			
			,chargedMultiplicity_(-9999.)
			,neutralMultiplicity_(-9999.)
			,muonMultiplicity_(-9999.)			
			,HFHadronMultiplicity_(-9999.)
			,HFEMMultiplicity_(-9999.)
			{;}

		~CatPFJet() {;}

		Float_t chargedHadronEnergyFraction() const { return chargedHadronEnergyFraction_; }
		Float_t neutralHadronEnergyFraction() const { return neutralHadronEnergyFraction_; }
		Float_t chargedEmEnergyFraction() const { return chargedEmEnergyFraction_; }
		Float_t chargedMuEnergyFraction() const { return chargedMuEnergyFraction_; }
		Float_t neutralEmEnergyFraction() const { return neutralEmEnergyFraction_; }		
		Float_t HFHadronEnergyFraction() const { return HFHadronEnergyFraction_; }
		Float_t HFEMEnergyFraction() const { return HFEMEnergyFraction_; }		
		Float_t chargedMultiplicity() const { return chargedMultiplicity_; }
		Float_t neutralMultiplicity() const { return neutralMultiplicity_; }
		Float_t muonMultiplicity() const { return muonMultiplicity_; }	  
		Float_t HFHadronMultiplicity() const { return HFHadronMultiplicity_; }
		Float_t HFEMMultiplicity() const { return HFEMMultiplicity_; }

		virtual TString typeName() const { return "CatPFJet"; }

		void setChargedHadronEnergyFraction(Float_t chargedHadronEnergyFraction) { chargedHadronEnergyFraction_ = chargedHadronEnergyFraction; }
		void setNeutralHadronEnergyFraction(Float_t neutralHadronEnergyFraction) { neutralHadronEnergyFraction_ = neutralHadronEnergyFraction; }
		void setChargedEmEnergyFraction(Float_t chargedEmEnergyFraction) { chargedEmEnergyFraction_ = chargedEmEnergyFraction; }
		void setChargedMuEnergyFraction(Float_t chargedMuEnergyFraction) { chargedMuEnergyFraction_ = chargedMuEnergyFraction; }
		void setNeutralEmEnergyFraction(Float_t neutralEmEnergyFraction) { neutralEmEnergyFraction_ = neutralEmEnergyFraction; }		
		void setHFHadronEnergyFraction(Float_t HFHadronEnergyFraction) { HFHadronEnergyFraction_ = HFHadronEnergyFraction; }
		void setHFEMEnergyFraction(Float_t HFEMEnergyFraction) { HFEMEnergyFraction_ = HFEMEnergyFraction; }		
		void setChargedMultiplicity(Float_t chargedMultiplicity) { chargedMultiplicity_ = chargedMultiplicity; }
		void setNeutralMultiplicity(Float_t neutralMultiplicity) { neutralMultiplicity_ = neutralMultiplicity; }
		void setMuonMultiplicity(Float_t muonMultiplicity) { muonMultiplicity_ = muonMultiplicity; }		
		void setHFHadronMultiplicity(Float_t HFHadronMultiplicity) { HFHadronMultiplicity_ = HFHadronMultiplicity; }
		void setHFEMMultiplicity(Float_t HFEMMultiplicity) { HFEMMultiplicity_ = HFEMMultiplicity; }

		friend std::ostream& operator<< (std::ostream& stream, const CatPFJet& jet)
		{
			stream << "CatPFJet - Charge=" << setw(2) << jet.charge() << " (Et,eta,phi)=("<< setw(10) << jet.Et() <<","<< setw(10) << jet.Eta() <<","<< setw(10) << jet.Phi() << ")"
				<< " vertex(x,y,z)=("<< jet.vx() <<","<< jet.vy() <<","<< jet.vz() << ")";
			return stream;
		};


	private:

		Float_t chargedHadronEnergyFraction_;
		Float_t neutralHadronEnergyFraction_;
		Float_t chargedEmEnergyFraction_;
		Float_t chargedMuEnergyFraction_;
		Float_t neutralEmEnergyFraction_;
		Float_t HFHadronEnergyFraction_;
		Float_t HFEMEnergyFraction_;
		Float_t chargedMultiplicity_;
		Float_t neutralMultiplicity_;
		Float_t muonMultiplicity_;
		Float_t HFHadronMultiplicity_;
		Float_t HFEMMultiplicity_;

		ClassDef (CatPFJet,2);
	};
}

#endif
