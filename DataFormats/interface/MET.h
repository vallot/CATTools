#ifndef CATTools_DataFormats_MET_h
#define CATTools_DataFormats_MET_h

#include "../interface/Particle.h"

#include "Rtypes.h"
#include "TObject.h"

using namespace std;

namespace cat
{
	class MET : public Particle
	{
	public:

		MET() :
			Particle()
      ,METType_(0)
			,sumEt_(-9999.)
			{;}

		MET(const MET& met) :
			Particle(met)
			,METType_(met.METType_)
			,sumEt_(met.sumEt_)
			{;}

		MET(Double_t px, Double_t py, Double_t pz, Double_t e) :
			Particle(px,py,pz,e)
			,METType_(0)
			,sumEt_(-9999.)
			{;}
	
		MET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			Particle(px,py,pz,e,vtx_x,vtx_y,vtx_z)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		MET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
			Particle(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		MET(const TLorentzVector &momentum) :
			Particle(momentum)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		MET(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
			Particle(momentum, vertex, type, charge)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		~MET() {;}
		
		Int_t METType() const { return METType_; }
		Float_t sumEt() const { return sumEt_; }
		virtual TString typeName() const { return "MET"; }

		void setMETType(Int_t METType) { METType_ = METType; }
		void setSumEt(Float_t sumEt) { sumEt_ = sumEt; }

		friend std::ostream& operator<< (std::ostream& stream, const MET& met)
		{
			stream << "MET  (Pt,Px,Py)=("<< met.Pt() <<","<< met.Px() <<","<< met.Py() << ")";
			return stream;
		};


	private:

		Int_t METType_;
		Float_t sumEt_;

		ClassDef (MET,2);
	};
}

#endif
