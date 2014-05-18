#ifndef CatLepton_h
#define CatLepton_h

#include <map>
#include <iostream>
#include <TLorentzVector.h>


#include "../interface/CatParticle.h"

using namespace std;

namespace cat
{
        class CatLepton : public CatParticle
	{
    
	public:
	        CatLepton() :
			CatParticle(),
			ip3d_(-9999.),
			ip3dErr_(-9999.),
			d0_(-9999.),
			d0Error_(-9999.),
			dz_(-9999.),
			dzError_(-9999.),
			//    dB_(-9999.),
			//    dBError_(-9999.),
			trackIso03_(-9999.),
			ecalIso03_(-9999.),
			hcalIso03_(-9999.),
			trackIso04_(-9999.),
			ecalIso04_(-9999.),
			hcalIso04_(-9999.),
			chargedHadronIso03_(-9999),
			puChargedHadronIso03_(-9999.),
			photonIso03_(-9999),
			neutralHadronIso03_(-9999),
			chargedHadronIso04_(-9999),
			puChargedHadronIso04_(-9999.),
			photonIso04_(-9999),
			neutralHadronIso04_(-9999)
			{;}
                
                CatLepton(const CatLepton& lepton) :
                        CatParticle(lepton),
                        ip3d_(lepton.ip3d_),
                        ip3dErr_(lepton.ip3dErr_),
                        d0_(lepton.d0_),
                        d0Error_(lepton.d0Error_),
                        dz_(lepton.dz_),
                        dzError_(lepton.dzError_),
                        //    dB_(lepton.dB_),
                        //    dBError_(lepton.dBError_),
                        trackIso03_(lepton.trackIso03_),
                        ecalIso03_(lepton.ecalIso03_),
                        hcalIso03_(lepton.hcalIso03_),
                        trackIso04_(lepton.trackIso04_),
                        ecalIso04_(lepton.ecalIso04_),
                        hcalIso04_(lepton.hcalIso04_),
                        chargedHadronIso03_(lepton.chargedHadronIso03_),
                        puChargedHadronIso03_(lepton.puChargedHadronIso03_),
                        photonIso03_(lepton.photonIso03_),
                        neutralHadronIso03_(lepton.neutralHadronIso03_),
                        chargedHadronIso04_(lepton.chargedHadronIso04_),
                        puChargedHadronIso04_(lepton.puChargedHadronIso04_),
                        photonIso04_(lepton.photonIso04_),
                        neutralHadronIso04_(lepton.neutralHadronIso04_)
                        {;}    

                CatLepton(Double_t px, Double_t py, Double_t pz, Double_t e) :
                        CatParticle(px,py,pz,e),
                        ip3d_(-9999.),
                        ip3dErr_(-9999.),
                        d0_(-9999.),
                        d0Error_(-9999.),
                        dz_(-9999.),
                        dzError_(-9999.),
                        //    dB_(-9999.),
                        //    dBError_(-9999.),
                        trackIso03_(-9999.),
                        ecalIso03_(-9999.),
                        hcalIso03_(-9999.),
                        trackIso04_(-9999.),
                        ecalIso04_(-9999.),
                        hcalIso04_(-9999.),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
                        {;}

                CatLepton(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
                        CatParticle(px,py,pz,e,vtx_x,vtx_y,vtx_z),
                        ip3d_(-9999.),
                        ip3dErr_(-9999.),
                        d0_(-9999.),
                        d0Error_(-9999.),
                        dz_(-9999.),
                        dzError_(-9999.),
                        //    dB_(-9999.),
                        //    dBError_(-9999.),
                        trackIso03_(-9999.),
                        ecalIso03_(-9999.),
                        hcalIso03_(-9999.),
                        trackIso04_(-9999.),
                        ecalIso04_(-9999.),
                        hcalIso04_(-9999.),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
                        {;}

                CatLepton(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Int_t charge) :
                        CatParticle(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge),
                        ip3d_(-9999.),
                        ip3dErr_(-9999.),
                        d0_(-9999.),
                        d0Error_(-9999.),
                        dz_(-9999.),
                        dzError_(-9999.),
                        //    dB_(-9999.),
                        //    dBError_(-9999.),
                        trackIso03_(-9999.),
                        ecalIso03_(-9999.),
                        hcalIso03_(-9999.),
                        trackIso04_(-9999.),
                        ecalIso04_(-9999.),
                        hcalIso04_(-9999.),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
                        {;}

                CatLepton(const TLorentzVector &momentum) :
                        CatParticle(momentum),
                        ip3d_(-9999.),
                        ip3dErr_(-9999.),
                        d0_(-9999.),
                        d0Error_(-9999.),
                        dz_(-9999.),
                        dzError_(-9999.),
                        //    dB_(-9999.),
                        //    dBError_(-9999.),
                        trackIso03_(-9999.),
                        ecalIso03_(-9999.),
                        hcalIso03_(-9999.),
                        trackIso04_(-9999.),
                        ecalIso04_(-9999.),
                        hcalIso04_(-9999.),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
                        {;}

                CatLepton(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
                        CatParticle(momentum,vertex,type,charge),
                        ip3d_(-9999.),
                        ip3dErr_(-9999.),
                        d0_(-9999.),
                        d0Error_(-9999.),
                        dz_(-9999.),
                        dzError_(-9999.),
                        //    dB_(-9999.),
                        //    dBError_(-9999.),
                        trackIso03_(-9999.),
                        ecalIso03_(-9999.),
                        hcalIso03_(-9999.),
                        trackIso04_(-9999.),
                        ecalIso04_(-9999.),
                        hcalIso04_(-9999.),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
                        {;}
 
	        ~CatLepton() {;}
    
    
	public:
		virtual TString typeName() const { return "CatLepton"; }
    
                Float_t ip3d() const { return ip3d_; }
                Float_t ip3dError() const { return ip3dErr_; }
		Float_t d0() const { return d0_; }
		Float_t d0Error()const { return d0Error_; }
		Float_t dz()const { return dz_; }
		Float_t dzError()const { return dzError_; }
    //		Float_t dB() const { return dB_; }
    //		Float_t dBError() const { return dBError_; }
    
                //detector based isolation
		Float_t ecalIso(unsigned int cone=3) const
		{
			if(cone == 3) return ecalIso03_;
			else if(cone == 4) return ecalIso04_;
			else cout<<"Bad Cone Size! It returns 9999."<<endl;
			return 9999.;
		}
		Float_t hcalIso(unsigned int cone=3) const
		{
			if(cone == 3) return hcalIso03_;
			else if(cone == 4) return hcalIso04_;
			else cout<<"Bad Cone Size! It returns 9999."<<endl;
			return 9999.;
  	        }
		Float_t caloIso(unsigned int cone=3) const
		{
			return (ecalIso(cone) + hcalIso(cone));
		}
		Float_t trackIso(unsigned int cone=3) const
		{
			if(cone == 3) return trackIso03_;
			else if(cone == 4) return trackIso04_;
			else cout<<"Bad Cone Size! It returns 9999."<<endl;
			return 9999.;
		}
		Float_t absDetIso(unsigned int cone=3) const
		{
			return (trackIso(cone) + hcalIso(cone) + ecalIso(cone));
		}
                Float_t relDetIso(unsigned int cone=3) const
                {
                  double relIso = absDetIso(cone)/((TLorentzVector)(*this)).Pt();
                  return relIso;
                }
    
                //particle based isolation
		Float_t chargedHadronIso(unsigned int cone=3) const 
                { 
                  if(cone == 3) return chargedHadronIso03_;
                  else if(cone == 4) return chargedHadronIso04_;
                  else cout <<"Bad Cone Size! It returns 9999."<<endl; 
                  return 9999; 
                }
                Float_t puChargedHadronIso(unsigned int cone=3) const
                { 
                  if(cone == 3) return puChargedHadronIso03_;
                  else if(cone == 4) return puChargedHadronIso04_;
                  else cout <<"Bad Cone Size! It returns 9999."<<endl;
                  return 9999;              
                }
                Float_t photonIso(unsigned int cone=3) const
                { 
                  if(cone == 3) return photonIso03_;
                  else if(cone == 4) return photonIso04_;
                  else cout <<"Bad Cone Size! It returns 9999."<<endl;
                  return 9999;              
                }
                Float_t neutralHadronIso(unsigned int cone=3) const
                {
                  if(cone == 3) return neutralHadronIso03_;
                  else if(cone == 4) return neutralHadronIso04_;
                  else cout <<"Bad Cone Size! It returns 9999."<<endl;
                  return 9999;              
                }
                Float_t absPfIso(unsigned int cone=3, float dBetaFactor=0) const
                {
                  // in this case, dbeta correction is asked, but 
                  // the input for this correction is not available. 
                  // better returning an unphysical result than applying a wrong correction.
                  if(dBetaFactor>0 && puChargedHadronIso(cone)<0) return -1;

                  double neutralIso = neutralHadronIso(cone) + photonIso(cone);
                  double corNeutralIso = neutralIso - dBetaFactor * puChargedHadronIso(cone);
                  double charged = chargedHadronIso(cone);

                  return charged + ( corNeutralIso>0 ? corNeutralIso : 0 ) ;
                }
		Float_t relPfIso(unsigned int cone=3, float dBetaFactor=0) const 
                { 
                  double relIso = absPfIso(cone, dBetaFactor)/((TLorentzVector)(*this)).Pt();
                  return relIso;
                }
    
		//setters
                void setIp3d(Float_t x) { ip3d_ = x; }
                void setIp3dError(Float_t x) { ip3dErr_ = x; }
		void setD0(Float_t x) { d0_ = x; }
		void setD0Error(Float_t d0Error) { d0Error_ = d0Error; }
		void setDz(Float_t x) { dz_ = x; }
		void setDzError(Float_t x) { dzError_ = x; }
		//void setDB(Float_t dB) { dB_ = dB; }
		//void setDBError(Float_t dBError) { dBError_ = dBError; }
   
                //detector based isolation
                void setIsoR03_trackIso(Float_t isoR03_trackIso) { trackIso03_ = isoR03_trackIso; }
                void setIsoR03_ecalIso(Float_t isoR03_ecalIso)   { ecalIso03_ = isoR03_ecalIso; }
                void setIsoR03_hcalIso(Float_t isoR03_hcalIso)   { hcalIso03_ = isoR03_hcalIso; }
                void setIsoR04_trackIso(Float_t isoR04_trackIso) { trackIso04_ = isoR04_trackIso; }
                void setIsoR04_ecalIso(Float_t isoR04_ecalIso)   { ecalIso04_ = isoR04_ecalIso; }
                void setIsoR04_hcalIso(Float_t isoR04_hcalIso)   { hcalIso04_ = isoR04_hcalIso; }
 
                //particle based isolation
                void setIsoR03_ChargedHadronIso(Float_t chargedHadronIso){ chargedHadronIso03_ = chargedHadronIso; }
                void setIsoR03_PuChargedHadronIso(Float_t iso) { puChargedHadronIso03_ = iso; }
                void setIsoR03_PhotonIso(Float_t photonIso){ photonIso03_ = photonIso; }
                void setIsoR03_NeutralHadronIso(Float_t neutralHadronIso){ neutralHadronIso03_ = neutralHadronIso; }
                void setIsoR04_ChargedHadronIso(Float_t chargedHadronIso){ chargedHadronIso04_ = chargedHadronIso; }
                void setIsoR04_PuChargedHadronIso(Float_t iso) { puChargedHadronIso04_ = iso; }
                void setIsoR04_PhotonIso(Float_t photonIso){ photonIso04_ = photonIso; }
                void setIsoR04_NeutralHadronIso(Float_t neutralHadronIso){ neutralHadronIso04_ = neutralHadronIso; }
	 
        protected:
      
                //TrackProperties=====================================
	        Float_t ip3d_;                             // 3D impact parameter
	        Float_t ip3dErr_;                          // error on ip3d_
	        Float_t d0_;                         	     // transverse impact parameter (wrt to PV)
	        Float_t d0Error_;                          // error on d0_
	        Float_t dz_;                               // longitudinal impact parameter (wrt to PV)
	        Float_t dzError_;                          // error on dz_ 
	      
	        // If this was not the case, dB is calculated wrt the beamspot and edb = -1 all the time
	        //Float_t dB_;                             // dB from PAT muon
	        //Float_t dBError_;                        // dBError from PAT muon
	      
	        //Isolation ============================================
                //detector based isolation
	        Float_t trackIso03_;                        // track iso deposit with electron footprint removed
	        Float_t ecalIso03_;                // ecal iso deposit with electron footprint removed
	        Float_t hcalIso03_;           // hcal depht 1 iso deposit with electron footprint removed
	        Float_t trackIso04_;                        // track iso deposit with electron footprint removed
	        Float_t ecalIso04_;                // ecal iso deposit with electron footprint removed
	        Float_t hcalIso04_;           // hcal depht 1 iso deposit with electron footprint removed

                //particle based isolation
                Float_t chargedHadronIso03_;                 // isolation calculated with only the charged hadron candidates
	        Float_t puChargedHadronIso03_;               // isolation calculated with only the pile-up charged hadron candidates
	        Float_t photonIso03_;                        // isolation calculated with only the gamma candidates
	        Float_t neutralHadronIso03_;	               // isolation calculated with only the neutral hadron candidates
	        Float_t chargedHadronIso04_;                 // isolation calculated with only the charged hadron candidates
	        Float_t puChargedHadronIso04_;               // isolation calculated with only the pile-up charged hadron candidates
                Float_t photonIso04_;                        // isolation calculated with only the gamma candidates
                Float_t neutralHadronIso04_;                     // isolation calculated with only the neutral hadron candidates
      
                ClassDef (CatLepton,3);
        };
}

#endif
