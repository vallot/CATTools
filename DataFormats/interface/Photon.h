#ifndef CATTools_DataFormats_Photon_h
#define CATTools_DataFormats_Photon_h

#include <map>
#include <iostream>
#include <TLorentzVector.h>

#include "../interface/Particle.h"

using namespace std;

namespace cat
{
	class Photon : public Particle 
	{

	public:
		Photon() :
			Particle(),
                        sigmaIetaIeta_(-9999.),
                        hadronicOverEm_(-9999.),
                        hasPixelSeed_(-9999.),
                        passelectronveto_(false),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
			{;}
	
		Photon(const Photon& photon) :
			Particle(photon),
                        sigmaIetaIeta_(photon.sigmaIetaIeta_),
                        hadronicOverEm_(photon.hadronicOverEm_),
                        hasPixelSeed_(photon.hasPixelSeed_),
                        passelectronveto_(photon.passelectronveto_),
                        chargedHadronIso03_(photon.chargedHadronIso03_),
                        puChargedHadronIso03_(photon.puChargedHadronIso03_),
                        photonIso03_(photon.photonIso03_),
                        neutralHadronIso03_(photon.neutralHadronIso03_),
                        chargedHadronIso04_(photon.chargedHadronIso04_),
                        puChargedHadronIso04_(photon.puChargedHadronIso04_),
                        photonIso04_(photon.photonIso04_),
                        neutralHadronIso04_(photon.neutralHadronIso04_)
			{;}
                Photon(Double_t px, Double_t py, Double_t pz, Double_t e) :
                        Particle(px,py,pz,e),
                        sigmaIetaIeta_(-9999.),
                        hadronicOverEm_(-9999.),
                        hasPixelSeed_(-9999.),
                        passelectronveto_(false),
                        chargedHadronIso03_(-9999),
                        puChargedHadronIso03_(-9999.),
                        photonIso03_(-9999),
                        neutralHadronIso03_(-9999),
                        chargedHadronIso04_(-9999),
                        puChargedHadronIso04_(-9999.),
                        photonIso04_(-9999),
                        neutralHadronIso04_(-9999)
                        {;}
      	
		~Photon() {;}
  
                Float_t sigmaIetaIeta() const { return sigmaIetaIeta_; } 
                Float_t hadronicOverEm() const { return hadronicOverEm_; }
                Float_t hasPixelSeed() const { return hasPixelSeed_; }
                Bool_t passelectronveto() const { return passelectronveto_; }

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
 
                void setSigmaIetaIeta( Float_t sigmaIetaIeta ) { sigmaIetaIeta_ = sigmaIetaIeta; }               	
                void setHadronicOverEm( Float_t hadronicOverEm ) { hadronicOverEm_ = hadronicOverEm; }
                void setHasPixelSeed( Float_t hasPixelSeed ) { hasPixelSeed_ = hasPixelSeed; }
                void setPasselectronveto( Bool_t passelectronveto ) { passelectronveto_ = passelectronveto; }

                //particle based isolation
                void setIsoR03_ChargedHadronIso(Float_t chargedHadronIso){ chargedHadronIso03_ = chargedHadronIso; }
                void setIsoR03_PuChargedHadronIso(Float_t iso) { puChargedHadronIso03_ = iso; }
                void setIsoR03_PhotonIso(Float_t photonIso){ photonIso03_ = photonIso; }
                void setIsoR03_NeutralHadronIso(Float_t neutralHadronIso){ neutralHadronIso03_ = neutralHadronIso; }
                void setIsoR04_ChargedHadronIso(Float_t chargedHadronIso){ chargedHadronIso04_ = chargedHadronIso; }
                void setIsoR04_PuChargedHadronIso(Float_t iso) { puChargedHadronIso04_ = iso; }
                void setIsoR04_PhotonIso(Float_t photonIso){ photonIso04_ = photonIso; }
                void setIsoR04_NeutralHadronIso(Float_t neutralHadronIso){ neutralHadronIso04_ = neutralHadronIso; }

	private:

                Float_t sigmaIetaIeta_;	
                Float_t hadronicOverEm_;
                Float_t hasPixelSeed_;	
                Bool_t passelectronveto_;

                //Isolation ============================================
                //particle based isolation
                Float_t chargedHadronIso03_;                 // isolation calculated with only the charged hadron candidates
                Float_t puChargedHadronIso03_;               // isolation calculated with only the pile-up charged hadron candidates
                Float_t photonIso03_;                        // isolation calculated with only the gamma candidates
                Float_t neutralHadronIso03_;                   // isolation calculated with only the neutral hadron candidates
                Float_t chargedHadronIso04_;                 // isolation calculated with only the charged hadron candidates
                Float_t puChargedHadronIso04_;               // isolation calculated with only the pile-up charged hadron candidates
                Float_t photonIso04_;                        // isolation calculated with only the gamma candidates
                Float_t neutralHadronIso04_;                     // isolation calculated with only the neutral hadron candidates

		ClassDef (Photon,2);
	};
}

#endif
