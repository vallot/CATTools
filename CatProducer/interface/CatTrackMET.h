#ifndef CatTrackMET_h
#define CatTrackMET_h

#include "../interface/CatParticle.h"
#include "../interface/CatMET.h"

#include "Rtypes.h"
#include "TObject.h"

using namespace std;

namespace cat
{
	class CatTrackMET : public CatMET
	{
	public:

		CatTrackMET() :
		  CatMET()
		
		  {;}

		CatTrackMET(const CatTrackMET& met) :
		  CatMET(met)
		  
		  {;}
		  

		  CatTrackMET(const CatMET& met) :
		    CatMET(met)
		   
		    {;}


		  CatTrackMET(Double_t px, Double_t py, Double_t pz, Double_t e) :
		    CatMET(px,py,pz,e)
		    
		    {;}
	
		CatTrackMET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
		  CatMET(px,py,pz,e,vtx_x,vtx_y,vtx_z)
		  		
		  {;}

		CatTrackMET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
		  CatMET(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
		 
		  {;}

		CatTrackMET(const TLorentzVector &momentum) :
		  CatMET(momentum)
		
		  {;}

		CatTrackMET(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
		  CatMET(momentum, vertex, type, charge)
		  		
		  {;}

		~CatTrackMET() {;}

		
	        virtual TString typeName() const { return "CatTrackMET"; }

		
		friend std::ostream& operator<< (std::ostream& stream, const CatTrackMET& met)
		{
			stream << "CatTrackMET  (Pt,Px,Py)=("<< met.Pt() <<","<< met.Px() <<","<< met.Py() << ")";
			return stream;
		};


	private:

		

		ClassDef (CatTrackMET,2);
	};
}

#endif
