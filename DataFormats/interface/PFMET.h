#ifndef CATTools_DataFormats_PFMET_h
#define CATTools_DataFormats_PFMET_h

#include "../interface/Particle.h"
#include "../interface/MET.h"

#include "Rtypes.h"
#include "TObject.h"

using namespace std;

namespace cat
{
	class PFMET : public MET
	{
	public:

		PFMET() :
		  MET()
		  ,NeutralEMFraction_(-9999.)
		  ,NeutralHadEtFraction_(-9999.)
		  ,ChargedEMEtFraction_(-9999.)
		  ,ChargedHadEtFraction_(-9999.)
		  ,MuonEtFraction_(-9999.)
		  ,Type6EtFraction_(-9999.)
		  ,Type7EtFraction_(-9999.)
		  {;}

		PFMET(const PFMET& met) :
		  MET(met)
		  ,NeutralEMFraction_(met.NeutralEMFraction_)
		  ,NeutralHadEtFraction_(met.NeutralHadEtFraction_)
		  ,ChargedEMEtFraction_(met.ChargedEMEtFraction_)
		  ,ChargedHadEtFraction_(met.ChargedHadEtFraction_)
		  ,MuonEtFraction_(met.MuonEtFraction_)
		  ,Type6EtFraction_(met.Type6EtFraction_)
		  ,Type7EtFraction_(met.Type7EtFraction_)
		  {;}
		  

		  PFMET(const MET& met) :
		    MET(met)
		    ,NeutralEMFraction_(-9999.)
		    ,NeutralHadEtFraction_(-9999.)
		    ,ChargedEMEtFraction_(-9999.)
		    ,ChargedHadEtFraction_(-9999.)
		    ,MuonEtFraction_(-9999.)
		    ,Type6EtFraction_(-9999.)
		    ,Type7EtFraction_(-9999.)
		    {;}


		  PFMET(Double_t px, Double_t py, Double_t pz, Double_t e) :
		    MET(px,py,pz,e)
		    ,NeutralEMFraction_(-9999.)
		    ,NeutralHadEtFraction_(-9999.)
		    ,ChargedEMEtFraction_(-9999.)
		    ,ChargedHadEtFraction_(-9999.)
		    ,MuonEtFraction_(-9999.)
		    ,Type6EtFraction_(-9999.)
		    ,Type7EtFraction_(-9999.)
		    {;}
	
		PFMET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
		  MET(px,py,pz,e,vtx_x,vtx_y,vtx_z)
		  ,NeutralEMFraction_(-9999.)
		  ,NeutralHadEtFraction_(-9999.)
		  ,ChargedEMEtFraction_(-9999.)
		  ,ChargedHadEtFraction_(-9999.)
		  ,MuonEtFraction_(-9999.)
		  ,Type6EtFraction_(-9999.)
		  ,Type7EtFraction_(-9999.)			
		  {;}

		PFMET(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
		  MET(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
		  ,NeutralEMFraction_(-9999.)
		  ,NeutralHadEtFraction_(-9999.)
		  ,ChargedEMEtFraction_(-9999.)
		  ,ChargedHadEtFraction_(-9999.)
		  ,MuonEtFraction_(-9999.)
		  ,Type6EtFraction_(-9999.)
		  ,Type7EtFraction_(-9999.)
		  {;}

		PFMET(const TLorentzVector &momentum) :
		  MET(momentum)
		  ,NeutralEMFraction_(-9999.)
		  ,NeutralHadEtFraction_(-9999.)
		  ,ChargedEMEtFraction_(-9999.)
		  ,ChargedHadEtFraction_(-9999.)
		  ,MuonEtFraction_(-9999.)
		  ,Type6EtFraction_(-9999.)
		  ,Type7EtFraction_(-9999.)
		  {;}

		PFMET(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
		  MET(momentum, vertex, type, charge)
		  ,NeutralEMFraction_(-9999.)
		  ,NeutralHadEtFraction_(-9999.)
		  ,ChargedEMEtFraction_(-9999.)
		  ,ChargedHadEtFraction_(-9999.)
		  ,MuonEtFraction_(-9999.)
		  ,Type6EtFraction_(-9999.)
		  ,Type7EtFraction_(-9999.)			
		  {;}

		~PFMET() {;}

		Double_t NeutralEMFraction() const { return NeutralEMFraction_; }
		Double_t NeutralHadEtFraction() const { return NeutralHadEtFraction_; }
		Double_t ChargedEMEtFraction() const { return ChargedEMEtFraction_; }
		Double_t ChargedHadEtFraction() const { return ChargedHadEtFraction_; }
		Double_t MuonEtFraction() const { return MuonEtFraction_; }
		Double_t Type6EtFraction() const { return Type6EtFraction_; }
		Double_t Type7EtFraction() const { return Type7EtFraction_; }

	
		//TObject* genMET() const { return genMET_.GetObject(); }
		virtual TString typeName() const { return "PFMET"; }

		void setPFMETFraction(
				      Double_t NeutralEMFraction
				      ,Double_t NeutralHadEtFraction
				      ,Double_t ChargedEMEtFraction
				      ,Double_t ChargedHadEtFraction
				      ,Double_t MuonEtFraction
				      ,Double_t Type6EtFraction
				      ,Double_t Type7EtFraction
				      )
		{
		  
		  NeutralEMFraction_ = NeutralEMFraction;
		  NeutralHadEtFraction_ = NeutralHadEtFraction;
		  ChargedEMEtFraction_ = ChargedEMEtFraction;
		  ChargedHadEtFraction_ = ChargedHadEtFraction;
		  MuonEtFraction_ = MuonEtFraction;
		  Type6EtFraction_ = Type6EtFraction;
		  Type7EtFraction_ = Type7EtFraction;

		}


		friend std::ostream& operator<< (std::ostream& stream, const PFMET& met)
		{
			stream << "PFMET  (Pt,Px,Py)=("<< met.Pt() <<","<< met.Px() <<","<< met.Py() << ")";
			return stream;
		};


	private:

		Double_t NeutralEMFraction_;
		Double_t NeutralHadEtFraction_;
		Double_t ChargedEMEtFraction_;
		Double_t ChargedHadEtFraction_;
		Double_t MuonEtFraction_;
		Double_t Type6EtFraction_;
		Double_t Type7EtFraction_;

		ClassDef (PFMET,2);
	};
}

#endif
