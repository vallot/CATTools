#ifndef CatJPTJet_h
#define CatJPTJet_h

#include "../interface/CatParticle.h"
#include "../interface/CatJet.h"

#include "Rtypes.h"
#include "TObject.h"

#include <iostream>
#include <iomanip>

using namespace std;

namespace cat
{
	class CatJPTJet : public CatJet
	{

	public:
		CatJPTJet() :
			CatJet()
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n90_(-9999)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
		  ,chargedHadronEnergyFraction_(-9999.)
			{;}
		
		CatJPTJet(const CatJPTJet& jet) :
			CatJet(jet)
			,etaetaMoment_(jet.etaetaMoment_)
			,phiphiMoment_(jet.phiphiMoment_)
			,ecalEnergyFraction_(jet.ecalEnergyFraction_)
			,hcalEnergyFraction_(jet.hcalEnergyFraction_)
			,maxEInEmTowers_(jet.maxEInEmTowers_)
			,maxEInHadTowers_(jet.maxEInHadTowers_)
			,towersArea_(jet.towersArea_)
			,n90_(jet.n90_)
			,n60_(jet.n60_)
			,fHPD_(jet.fHPD_)
			,fRBX_(jet.fRBX_)
			,n90Hits_(jet.n90Hits_)
			,nHCALTowers_(jet.nHCALTowers_)
			,nECALTowers_(jet.nECALTowers_)
			,chargedMultiplicity_(jet.chargedMultiplicity_)
		  ,chargedHadronEnergyFraction_(jet.chargedHadronEnergyFraction_)
			{;}
	
		CatJPTJet(const CatJet& jet) :
			CatJet(jet)
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n90_(-9999)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
		  ,chargedHadronEnergyFraction_(-9999.)
			{;}
	
		CatJPTJet(Double_t px, Double_t py, Double_t pz, Double_t e) :
			CatJet(px,py,px,e)
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n90_(-9999)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
		  ,chargedHadronEnergyFraction_(-9999.)
			{;}
	
		CatJPTJet(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			CatJet(px,py,pz,e,vtx_x,vtx_y,vtx_z)
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n90_(-9999)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
		  ,chargedHadronEnergyFraction_(-9999.)
			{;}
	
		CatJPTJet(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
			CatJet(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
		  ,chargedHadronEnergyFraction_(-9999.)
			{;}
	
		CatJPTJet(const TLorentzVector &momentum) :
			CatJet(momentum)
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n90_(-9999)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
		  ,chargedHadronEnergyFraction_(-9999.)
			{;}
	
		CatJPTJet(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
			CatJet(momentum, vertex, type, charge)
			,etaetaMoment_(-9999.)
			,phiphiMoment_(-9999.)
			,ecalEnergyFraction_(-9999.)
			,hcalEnergyFraction_(-9999.)
			,maxEInEmTowers_(-9999.)
   	  ,maxEInHadTowers_(-9999.)
			,towersArea_(-9999.)
			,n90_(-9999)
			,n60_(-9999)
			,fHPD_(-9999)
			,fRBX_(-9999.)
			,n90Hits_(-9999.)
			,nHCALTowers_(-9999)
			,nECALTowers_(-9999)
			,chargedMultiplicity_(-9999)
			,chargedHadronEnergyFraction_(-9999.)
			{;}
		
		~CatJPTJet() {;}
		
		Float_t etaetaMoment() const { return etaetaMoment_; }
		Float_t phiphiMoment() const { return phiphiMoment_; }
		Float_t ecalEnergyFraction() const { return ecalEnergyFraction_; }
		Float_t hcalEnergyFraction() const { return hcalEnergyFraction_; }
		Float_t maxEInEmTowers() const { return maxEInEmTowers_; }
		Float_t maxEInHadTowers() const { return maxEInHadTowers_;}
		Float_t towersArea() const { return towersArea_;} 
		Int_t n90() const { return n90_; }
		Int_t n60() const { return n60_; }
		Float_t fHPD() const { return fHPD_; }
		Float_t fRBX() const { return fRBX_; }
		Float_t n90Hits() const { return n90Hits_; }
		Int_t nHCALTowers() const { return nHCALTowers_; }
		Int_t nECALTowers() const { return nECALTowers_; }
		Int_t chargedMultiplicity() const { return chargedMultiplicity_; }
		Float_t chargedHadronEnergyFraction () const {return chargedHadronEnergyFraction_;}

		virtual TString typeName() const { return "CatJPTJet"; }
		
		void setetaetaMoment(Float_t etaetaMoment) { etaetaMoment_ = etaetaMoment; }
		void setphiphiMoment(Float_t phiphiMoment) { phiphiMoment_ = phiphiMoment; }
		void setEcalEnergyFraction(Float_t ecalEnergyFraction) { ecalEnergyFraction_ = ecalEnergyFraction; }
		void setHcalEnergyFraction(Float_t hcalEnergyFraction) { hcalEnergyFraction_ = hcalEnergyFraction; }
		void setMaxEInEmTowers(Float_t maxEInEmTowers) { maxEInEmTowers_ = maxEInEmTowers; }
		void setMaxEInHadTowers(Float_t maxEInHadTowers) { maxEInHadTowers_ = maxEInHadTowers; }
		void setTowersArea(Float_t towersArea) {towersArea_ = towersArea; }
		void setN90(Int_t n90) { n90_ = n90; }
		void setN60(Int_t n60) { n60_ = n60; }
		void setfHPD(Float_t fHPD) { fHPD_ = fHPD; }
		void setfRBX(Float_t fRBX) { fRBX_ = fRBX; }
		void setn90Hits(Float_t n90Hits) { n90Hits_ = n90Hits; }
		void setnHCALTowers(Int_t nHCALTowers) { nHCALTowers_ = nHCALTowers; }
		void setnECALTowers(Int_t nECALTowers) { nECALTowers_ = nECALTowers; }
		void setChargedMultiplicity(Int_t chargedMultiplicity) { chargedMultiplicity_ = chargedMultiplicity; }
		void setchargedHadronEnergyFraction (Float_t chargedHadronEnergyFraction) { chargedHadronEnergyFraction_ = chargedHadronEnergyFraction; }
		
		friend std::ostream& operator<< (std::ostream& stream, const CatJPTJet& jet) {
			stream << "CatJPTJet - Charge=" << setw(2) << jet.charge() << " (Et,eta,phi)=("<< setw(10) << jet.Et() <<","<< setw(10) << jet.Eta() <<","<< setw(10) << jet.Phi() << ")"
					<< " vertex(x,y,z)=("<< jet.vx() <<","<< jet.vy() <<","<< jet.vz() << ")";
			return stream;
		};


	private:

		Float_t etaetaMoment_;					// Added to CaloJet since they seem to be empty for PF
		Float_t phiphiMoment_;					// Added to CaloJet since they seem to be empty for PF
		Float_t ecalEnergyFraction_;        // ECAL Energy Fraction
		Float_t hcalEnergyFraction_;        // HCAL Energy Fraction
		Float_t maxEInEmTowers_;
		Float_t maxEInHadTowers_;
		Float_t towersArea_; 
		Int_t n90_;                         // Number of constituents of the jet carrying 90% of tje jet energy
		Int_t n60_;                         // Number of constituents of the jet carrying 60% of tje jet energy
		Float_t fHPD_;
		Float_t fRBX_;
		Float_t n90Hits_;
		Int_t nHCALTowers_;
		Int_t nECALTowers_;
		Int_t chargedMultiplicity_;

		Float_t chargedHadronEnergyFraction_;
		// neutralHadronEnergyFraction_, chargedEmEnergyFraction_, neutralEmEnergyFraction_ all not filled in PAT, so removed

		ClassDef (CatJPTJet,2);
	};
}

#endif
