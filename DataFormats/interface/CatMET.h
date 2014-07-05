#ifndef CatMET_h
#define CatMET_h

#include "../interface/CatParticle.h"

#include "Rtypes.h"
#include "TObject.h"

using namespace std;

namespace cat
{
	class CatMET : public CatParticle
	{
	public:

		CatMET() :
			CatParticle()
      ,METType_(0)
			,sumEt_(-9999.)
			{;}

		CatMET(const CatMET& met) :
			CatParticle(met)
			,METType_(met.METType_)
			,sumEt_(met.sumEt_)
			{;}

		CatMET(Double_t px, Double_t py, Double_t pz, Double_t e) :
			CatParticle(px,py,pz,e)
			,METType_(0)
			,sumEt_(-9999.)
			{;}
	
		CatMET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			CatParticle(px,py,pz,e,vtx_x,vtx_y,vtx_z)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		CatMET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
			CatParticle(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		CatMET(const TLorentzVector &momentum) :
			CatParticle(momentum)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		CatMET(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
			CatParticle(momentum, vertex, type, charge)
			,METType_(0)
			,sumEt_(-9999.)
			{;}

		~CatMET() {;}
		
		Int_t METType() const { return METType_; }
		Float_t sumEt() const { return sumEt_; }
		virtual TString typeName() const { return "CatMET"; }

		void setMETType(Int_t METType) { METType_ = METType; }
		void setSumEt(Float_t sumEt) { sumEt_ = sumEt; }

		friend std::ostream& operator<< (std::ostream& stream, const CatMET& met)
		{
			stream << "CatMET  (Pt,Px,Py)=("<< met.Pt() <<","<< met.Px() <<","<< met.Py() << ")";
			return stream;
		};


	private:

		Int_t METType_;
		Float_t sumEt_;

		ClassDef (CatMET,2);
	};
}

#endif
