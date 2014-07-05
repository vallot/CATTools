#ifndef Vertex_h
#define Vertex_h

#include <string>
#include <iostream>

#include "Rtypes.h"
#include "TVector3.h"
#include "TRef.h"

using namespace std;

namespace cat
{
	class Vertex : public TVector3
	{

	public:
		Vertex() :
			TVector3()
			,isValid_(false)
			,isFake_(true)
			,chi2_(-9999.)
			,ndof_(-9999.)
			,tracksSize_(-9999)
			,xError_(-9999.)
			,yError_(-9999.)
			,zError_(-9999.)
			{;}
	
		Vertex(const Vertex& vertex) :
			TVector3(vertex)
			,isValid_(vertex.isValid_)
			,isFake_(vertex.isFake_)
			,chi2_(vertex.chi2_)
			,ndof_(vertex.ndof_)
			,tracksSize_(vertex.tracksSize_)
			,xError_(vertex.xError_)
			,yError_(vertex.yError_)
			,zError_(vertex.zError_)
			{;}
	
		Vertex(Float_t x, Float_t y, Float_t z) :
			TVector3(x,y,z)
			,isValid_(false)
			,isFake_(true)
			,chi2_(-9999.)
			,ndof_(-9999.)
			,tracksSize_(-9999)
			,xError_(-9999.)
			,yError_(-9999.)
			,zError_(-9999.)
			{;}

		~Vertex() {;}
	
		Bool_t isValid() const { return isValid_; }
		Bool_t isFake() const { return isFake_; }
		Float_t chi2() const { return chi2_; }
		Float_t ndof() const { return ndof_; }
		Float_t normalizedChi2() const
		{
			Float_t normalizedChi2_ = -9999;
			if(ndof_ != 0) normalizedChi2_ = chi2_/ndof_;
			return normalizedChi2_;
		}
		Int_t tracksSize() const { return tracksSize_; }

		Float_t xError() const { return xError_; }
		Float_t yError() const { return yError_; }
		Float_t zError() const { return zError_; }

		void setIsValid(Bool_t isValid) { isValid_ = isValid; }
		void setIsFake(Bool_t isFake) { isFake_ = isFake; }
		void setChi2(Float_t chi2) { chi2_ = chi2; }
		void setNdof(Float_t ndof) { ndof_ = ndof; }
		void setTracksSize(Int_t tracksSize) { tracksSize_ = tracksSize; }

		void setXError(Float_t xError) { xError_ = xError; }
		void setYError(Float_t yError) { yError_ = yError; }
		void setZError(Float_t zError) { zError_ = zError; }

	private:
	
		Bool_t isValid_;
		Bool_t isFake_;
		Float_t chi2_;			//	Not divided by ndof
		Float_t ndof_;
		Int_t tracksSize_;
		
		Float_t xError_;
		Float_t yError_;
		Float_t zError_;

		ClassDef (Vertex,2);
	};
}

#endif
