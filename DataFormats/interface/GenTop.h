#ifndef CATTools_DataFormats_GenTop_h
#define CATTools_DataFormats_GenTop_h
#include "MCParticle.h"

#include <Math/VectorUtil.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

namespace cat
{
	class GenTop : public MCParticle
	{
		public:
		GenTop() :
			MCParticle()
			,isHadronic_(false)
			,isLeptonic_(false)
			,W_ (MCParticle())
			,bquark_ (MCParticle())
			,quark_ (MCParticle())
			,quarkBar_ (MCParticle())
			,lepton_ (MCParticle())
			,neutrino_ (MCParticle())
			,production_ (string(""))
			{};

		GenTop(const MCParticle& top) :
			MCParticle(top)
			,isHadronic_(false)
			,isLeptonic_(false)
			,W_ (MCParticle())
			,bquark_ (MCParticle())
			,quark_ (MCParticle())
			,quarkBar_ (MCParticle())
			,lepton_ (MCParticle())
			,neutrino_ (MCParticle())
			,production_ (string(""))
			{};

		//prod1: lepton prod2: neutrino   or prod1: quark   prod2: quarkBar
		GenTop(const bool isLeptonic, MCParticle& top, MCParticle& W, MCParticle& b, MCParticle& prod1, MCParticle& prod2, string production = string("") )
		{
			(*this) = GenTop(top);
			W_ = W;
			bquark_ = b;
			production_ = production;
			if(isLeptonic)
			{
				isLeptonic_ = isLeptonic;
				lepton_ = prod1;
				neutrino_ = prod2;
			}
			else
			{
				isHadronic_ = true;
				quark_ = prod1;
				quarkBar_ = prod2;
			}
		};

		GenTop(const GenTop& top):MCParticle(top){
			//(*this) = (MCParticle) top;
			isHadronic_ = top.isHadronic_;
			isLeptonic_ = top.isLeptonic_;
			W_ = top.W_;
			bquark_ = top.bquark_;
			quark_ = top.quark_;
			quarkBar_ = top.quarkBar_;
			lepton_ = top.lepton_;
			neutrino_ = top.neutrino_;
			production_ = top.production_;
		};
		~GenTop(){;};

		Bool_t isHadronic() const { return isHadronic_; }
		Bool_t isLeptonic() const { return isLeptonic_; }
		Bool_t isLeptonicMu() const { if(isLeptonic_ && abs(lepton().type())==13) return true ; else return false ;}
		Bool_t isLeptonicEl() const { if(isLeptonic_ && abs(lepton().type())==11) return true ; else return false ;}
		Bool_t isLeptonicTau() const { if(isLeptonic_ && abs(lepton().type())==15) return true ; else return false ;}

		Bool_t isHadronicWellSeparated(float deltaMin = 1.0) const;
		Float_t DeltaRMinHadronicTop() const;
 
		void Production() { cout<<production_<<endl; }
		Int_t From() const { return this->motherType(); }
		Bool_t FromGluino() const { return this->motherType()==1000021 ? true:false; }
		Bool_t FromStop() const { return (abs(this->motherType())==1000006 || abs(this->motherType())==2000006 ) ? true:false; }
		Bool_t FromSbottom() const { return (abs(this->motherType())==1000005 || abs(this->motherType())==2000005 ) ? true:false; }
 
		const MCParticle W() const {return W_;};
		const MCParticle bquark() const {return bquark_;};
		const MCParticle quark() const  {return quark_;};
		const MCParticle quarkBar() const {return quarkBar_;};
		const MCParticle lepton() const {return lepton_;};
		const MCParticle neutrino() const {return neutrino_;};
		friend std::ostream& operator<< (std::ostream& stream, const GenTop& top)
		{
			stream<<"-------------------------------------------"<<endl;
			stream <<"Top ";
			if(top.isHadronic()) stream << "hadronic ";
			if(top.isLeptonic()) stream << "leptonic ";
			if(top.isLeptonicMu()) stream <<"muon ";
			if(top.isLeptonicEl()) stream <<"electron ";
			if(top.isLeptonicTau()) stream <<"tau ";
			stream<<endl;
			stream<<top<<endl;
			stream<<"-------------------------------------------"<<endl;
			return stream ;
		};
 
	private:
 
		Bool_t isHadronic_;
		Bool_t isLeptonic_;
		MCParticle W_;
		MCParticle bquark_;
		MCParticle quark_;
		MCParticle quarkBar_;
		MCParticle lepton_;
		MCParticle neutrino_;
		string production_;

		ClassDef (GenTop,1);
	};
}

#endif
