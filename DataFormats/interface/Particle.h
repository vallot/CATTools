#ifndef Particle_h
#define Particle_h

#include <string>
#include <iostream>

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRef.h"

using namespace std;

namespace cat
{
	class Particle : public TLorentzVector
	{

	public:
	
		Particle() :
			TLorentzVector()
			,vertex_()
			,type_(0)
			,charge_(0.)
			,genParticleIndex_(-1)
			{;}

		Particle(const Particle& particle) :
			TLorentzVector(particle)
			,vertex_(particle.vertex_)
			,type_(particle.type_)
			,charge_(particle.charge_)
			,genParticleIndex_(particle.genParticleIndex_)
			{;}

		Particle(Double_t px, Double_t py, Double_t pz, Double_t e) :
			TLorentzVector(px,py,pz,e)
			,vertex_()
			,type_(0)
			,charge_(0.)
			,genParticleIndex_(-1)
			{;}

		Particle(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			TLorentzVector(px,py,pz,e)
			,vertex_(vtx_x,vtx_y,vtx_z)
			,type_(0)
			,charge_(0.)
			,genParticleIndex_(-1)
			{;}

		Particle(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z,Int_t type, Float_t charge) :
			TLorentzVector(px,py,pz,e)
			,vertex_(vtx_x,vtx_y,vtx_z)
			,type_(type)
			,charge_(charge)
			,genParticleIndex_(-1)
			{;}


		Particle(const TLorentzVector &momentum) :
			TLorentzVector(momentum)
			,vertex_()
			,type_(0)
			,charge_(0.)
			,genParticleIndex_(-1)
			{;}

		Particle(const TLorentzVector &momentum, const TVector3 &vertex) :
			TLorentzVector(momentum)
			,vertex_(vertex)
			,type_(0)
			,charge_(0.)
			,genParticleIndex_(-1)
			{;}

		Particle(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
			TLorentzVector(momentum)
			,vertex_(vertex)
			,type_(type)
			,charge_(charge)
			,genParticleIndex_(-1)
			{;}

		~Particle() {;}


		Double_t vx() const  { return vertex_.x(); }
		Double_t vy() const  { return vertex_.y(); }
		Double_t vz() const  { return vertex_.z(); }
		Int_t type() const  { return type_; }
		Float_t charge() const  { return charge_; }
		Int_t genParticleIndex() const { return genParticleIndex_; }
		virtual TString typeName() const { return "Particle"; }


		// FIXME setVx, setVy and setVz must modify the TLorentzVector ?
		void setVx(Double_t vx) { vertex_.SetX(vx); }
		void setVy(Double_t vy) { vertex_.SetY(vy); }
		void setVz(Double_t vz) { vertex_.SetZ(vz); }
		void setType(Int_t type) { type_ = type; }
		void setCharge(Int_t charge) { charge_ = charge; }
		void setGenParticleIndex(Int_t genParticleIndex) { genParticleIndex_ = genParticleIndex; }

		friend std::ostream& operator<< (std::ostream& stream, const Particle& part)
		{
			stream << "Type=" << part.type_ << "  Charge=" << part.charge_ << " (Et,eta,phi)=("<< part.Et() <<","<< part.Eta() <<","<< part.Phi() << ")"
				<< " vertex(x,y,z)=("<< part.vx() <<","<< part.vy() <<","<< part.vz() << ")";
			return stream;
		};


	protected:
	
		TVector3 vertex_;
		Int_t type_;
		Float_t charge_;
		Int_t genParticleIndex_;

		ClassDef (Particle,3);
	};
}

#endif
