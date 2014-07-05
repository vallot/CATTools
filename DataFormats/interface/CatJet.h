#ifndef CatJet_h
#define CatJet_h

#include "../interface/CatParticle.h"

#include "Rtypes.h"
#include "TObject.h"

#include <iostream>
#include <iomanip>
#include <map>

// Specific methods for PF and Calo can be found on:
// http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_5_2/doc/html/df/d60/DataFormats_2PatCandidates_2interface_2Jet_8h-source.html#l00212

using namespace std;

namespace cat
{
	class CatJet : public CatParticle
	{
	
	public:
		CatJet() :
			CatParticle()
			,jetType_(0)
			,nConstituents_(-9999)
			,jetArea_(-9999.)
			,maxDistance_(-9999.)
			,btag_jetBProbabilityBJetTags_(-9999.)
			,btag_jetProbabilityBJetTags_(-9999.)
			,btag_trackCountingHighPurBJetTags_(-9999.)
			,btag_trackCountingHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighPurBJetTags_(-9999.)
			,btag_combinedSecondaryVertexBJetTags_(-9999.)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(-9999.)
			,btag_combinedSecondaryVertexMVABJetTags_(-9999.)
			,btag_softMuonBJetTags_(-9999.)
			,btag_softMuonByPtBJetTags_(-9999.)
			,btag_softMuonByIP3dBJetTags_(-9999.)
			,btag_softElectronByPtBJetTags_(-9999.)
			,btag_softElectronByIP3dBJetTags_(-9999.)
			,btag_combinedCSVJPBJetTags_(-9999.)
			,btag_combinedCSVJPSLBJetTags_(-9999.)
			,btag_combinedCSVSLBJetTags_(-9999.)
			,btag_softPFElectronRetrainedBJetsTags_(-9999.)
			,btag_softPFMuonRetrainedBJetsTags_(-9999.)
			,partonFlavour_(-999)
			,isTopJet_(false)
			{;}
	
		CatJet(const CatJet& jet) :
			CatParticle(jet)
			,jetType_(jet.jetType_)
			,nConstituents_(jet.nConstituents_)
			,jetArea_(jet.jetArea_)
			,maxDistance_(jet.maxDistance_)
			,btag_jetBProbabilityBJetTags_(jet.btag_jetBProbabilityBJetTags_)
			,btag_jetProbabilityBJetTags_(jet.btag_jetProbabilityBJetTags_)
			,btag_trackCountingHighPurBJetTags_(jet.btag_trackCountingHighPurBJetTags_)
			,btag_trackCountingHighEffBJetTags_(jet.btag_trackCountingHighEffBJetTags_)
			,btag_simpleSecondaryVertexHighEffBJetTags_(jet.btag_simpleSecondaryVertexHighEffBJetTags_)
			,btag_simpleSecondaryVertexHighPurBJetTags_(jet.btag_simpleSecondaryVertexHighPurBJetTags_)
			,btag_combinedSecondaryVertexBJetTags_(jet.btag_combinedSecondaryVertexBJetTags_)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(jet.btag_combinedSecondaryVertexRetrainedBJetTags_)
			,btag_combinedSecondaryVertexMVABJetTags_(jet.btag_combinedSecondaryVertexMVABJetTags_)
			,btag_softMuonBJetTags_(jet.btag_softMuonBJetTags_)
			,btag_softMuonByPtBJetTags_(jet.btag_softMuonByPtBJetTags_)
			,btag_softMuonByIP3dBJetTags_(jet.btag_softMuonByIP3dBJetTags_)
			,btag_softElectronByPtBJetTags_(jet.btag_softElectronByPtBJetTags_)
			,btag_softElectronByIP3dBJetTags_(jet.btag_softElectronByIP3dBJetTags_)
			,btag_combinedCSVJPBJetTags_(jet.btag_combinedCSVJPBJetTags_)
			,btag_combinedCSVJPSLBJetTags_(jet.btag_combinedCSVJPSLBJetTags_)
			,btag_combinedCSVSLBJetTags_(jet.btag_combinedCSVSLBJetTags_)
			,btag_softPFElectronRetrainedBJetsTags_(jet.btag_softPFElectronRetrainedBJetsTags_)
			,btag_softPFMuonRetrainedBJetsTags_(jet.btag_softPFMuonRetrainedBJetsTags_)
			,partonFlavour_(jet.partonFlavour_)
			,isTopJet_(jet.isTopJet_)
		  {
			unsigned int size = sizeof(JetCorrName_)/sizeof(JetCorrName_[0]);
			for (unsigned int i=0; i<size; i++)
			{
				JetCorrName_[i] = jet.JetCorrName_[i];
			    	JetCorrValue_[i] = jet.JetCorrValue_[i];
			}
			for(std::map<std::string,float>::const_iterator it = jet.mistag_SF_.begin(); it != jet.mistag_SF_.end(); it++) {
				mistag_SF_[it->first] = it->second;
			}
			for(std::map<std::string,float>::const_iterator it = jet.btag_SF_.begin(); it != jet.btag_SF_.end(); it++) {
				btag_SF_[it->first] = it->second;
			}
			for(std::map<std::string,float>::const_iterator it = jet.mistag_SFerr_.begin(); it != jet.mistag_SFerr_.end(); it++) {
				mistag_SFerr_[it->first] = it->second;
			}
			for(std::map<std::string,float>::const_iterator it = jet.btag_SFerr_.begin(); it != jet.btag_SFerr_.end(); it++) {
				btag_SFerr_[it->first] = it->second;
			}
		}

		CatJet(Double_t px, Double_t py, Double_t pz, Double_t e) :
			CatParticle(px,py,pz,e)
			,jetType_(0)
			,nConstituents_(-9999)
			,jetArea_(-9999.)
			,maxDistance_(-9999.)
			,btag_jetBProbabilityBJetTags_(-9999.)
			,btag_jetProbabilityBJetTags_(-9999.)
			,btag_trackCountingHighPurBJetTags_(-9999.)
			,btag_trackCountingHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighPurBJetTags_(-9999.)
			,btag_combinedSecondaryVertexBJetTags_(-9999.)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(-9999.)
			,btag_combinedSecondaryVertexMVABJetTags_(-9999.)
			,btag_softMuonBJetTags_(-9999.)
			,btag_softMuonByPtBJetTags_(-9999.)
			,btag_softMuonByIP3dBJetTags_(-9999.)
			,btag_softElectronByPtBJetTags_(-9999.)
			,btag_softElectronByIP3dBJetTags_(-9999.)
			,btag_combinedCSVJPBJetTags_(-9999.)
			,btag_combinedCSVJPSLBJetTags_(-9999.)
			,btag_combinedCSVSLBJetTags_(-9999.)
			,btag_softPFElectronRetrainedBJetsTags_(-9999.)
			,btag_softPFMuonRetrainedBJetsTags_(-9999.)
			,partonFlavour_(-9999)
			,isTopJet_(false)
			{;}
	
		CatJet(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			CatParticle(px,py,pz,e,vtx_x,vtx_y,vtx_z)
			,jetType_(0)
			,nConstituents_(-9999)
			,jetArea_(-9999.)
			,maxDistance_(-9999.)
			,btag_jetBProbabilityBJetTags_(-9999.)
			,btag_jetProbabilityBJetTags_(-9999.)
			,btag_trackCountingHighPurBJetTags_(-9999.)
			,btag_trackCountingHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighPurBJetTags_(-9999.)
			,btag_combinedSecondaryVertexBJetTags_(-9999.)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(-9999.)
			,btag_combinedSecondaryVertexMVABJetTags_(-9999.)
			,btag_softMuonBJetTags_(-9999.)
			,btag_softMuonByPtBJetTags_(-9999.)
			,btag_softMuonByIP3dBJetTags_(-9999.)
			,btag_softElectronByPtBJetTags_(-9999.)
			,btag_softElectronByIP3dBJetTags_(-9999.)
			,btag_combinedCSVJPBJetTags_(-9999.)
			,btag_combinedCSVJPSLBJetTags_(-9999.)
			,btag_combinedCSVSLBJetTags_(-9999.)
			,btag_softPFElectronRetrainedBJetsTags_(-9999.)
			,btag_softPFMuonRetrainedBJetsTags_(-9999.)
			,partonFlavour_(-9999)
			,isTopJet_(false)
			{;}

		CatJet(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Float_t charge) :
			CatParticle(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
			,jetType_(0)
			,nConstituents_(-9999)
			,jetArea_(-9999.)
			,maxDistance_(-9999.)
			,btag_jetBProbabilityBJetTags_(-9999.)
			,btag_jetProbabilityBJetTags_(-9999.)
			,btag_trackCountingHighPurBJetTags_(-9999.)
			,btag_trackCountingHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighPurBJetTags_(-9999.)
			,btag_combinedSecondaryVertexBJetTags_(-9999.)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(-9999.)
			,btag_combinedSecondaryVertexMVABJetTags_(-9999.)
			,btag_softMuonBJetTags_(-9999.)
			,btag_softMuonByPtBJetTags_(-9999.)
			,btag_softMuonByIP3dBJetTags_(-9999.)
			,btag_softElectronByPtBJetTags_(-9999.)
			,btag_softElectronByIP3dBJetTags_(-9999.)
			,btag_combinedCSVJPBJetTags_(-9999.)
			,btag_combinedCSVJPSLBJetTags_(-9999.)
			,btag_combinedCSVSLBJetTags_(-9999.)
			,btag_softPFElectronRetrainedBJetsTags_(-9999.)
			,btag_softPFMuonRetrainedBJetsTags_(-9999.)
			,partonFlavour_(-9999)
			,isTopJet_(false)
			{;}

		CatJet(const TLorentzVector &momentum) :
			CatParticle(momentum)
			,jetType_(0)
			,nConstituents_(-9999)
			,jetArea_(-9999.)
			,maxDistance_(-9999.)
			,btag_jetBProbabilityBJetTags_(-9999.)
			,btag_jetProbabilityBJetTags_(-9999.)
			,btag_trackCountingHighPurBJetTags_(-9999.)
			,btag_trackCountingHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighPurBJetTags_(-9999.)
			,btag_combinedSecondaryVertexBJetTags_(-9999.)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(-9999.)
			,btag_combinedSecondaryVertexMVABJetTags_(-9999.)
			,btag_softMuonBJetTags_(-9999.)
			,btag_softMuonByPtBJetTags_(-9999.)
			,btag_softMuonByIP3dBJetTags_(-9999.)
			,btag_softElectronByPtBJetTags_(-9999.)
			,btag_softElectronByIP3dBJetTags_(-9999.)
			,btag_combinedCSVJPBJetTags_(-9999.)
			,btag_combinedCSVJPSLBJetTags_(-9999.)
			,btag_combinedCSVSLBJetTags_(-9999.)
			,btag_softPFElectronRetrainedBJetsTags_(-9999.)
			,btag_softPFMuonRetrainedBJetsTags_(-9999.)
			,partonFlavour_(-9999)
			,isTopJet_(false)
			{;}

		CatJet(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Float_t charge) :
			CatParticle(momentum, vertex, type, charge)
			,jetType_(0)
			,nConstituents_(-9999)
			,jetArea_(-9999.)
			,maxDistance_(-9999.)
			,btag_jetBProbabilityBJetTags_(-9999.)
			,btag_jetProbabilityBJetTags_(-9999.)
			,btag_trackCountingHighPurBJetTags_(-9999.)
			,btag_trackCountingHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighEffBJetTags_(-9999.)
			,btag_simpleSecondaryVertexHighPurBJetTags_(-9999.)
			,btag_combinedSecondaryVertexBJetTags_(-9999.)
			,btag_combinedSecondaryVertexRetrainedBJetTags_(-9999.)
			,btag_combinedSecondaryVertexMVABJetTags_(-9999.)
			,btag_softMuonBJetTags_(-9999.)
			,btag_softMuonByPtBJetTags_(-9999.)
			,btag_softMuonByIP3dBJetTags_(-9999.)
			,btag_softElectronByPtBJetTags_(-9999.)
			,btag_softElectronByIP3dBJetTags_(-9999.)
			,btag_combinedCSVJPBJetTags_(-9999.)
			,btag_combinedCSVJPSLBJetTags_(-9999.)
			,btag_combinedCSVSLBJetTags_(-9999.)
			,btag_softPFElectronRetrainedBJetsTags_(-9999.)
			,btag_softPFMuonRetrainedBJetsTags_(-9999.)
			,partonFlavour_(-9999)
			,isTopJet_(false)
			{;}

		~CatJet() {;}

		Int_t jetType() const { return jetType_; }
		Int_t nConstituents() const { return nConstituents_; }
		Float_t jetArea() const { return jetArea_; }
		Float_t maxDistance() const { return maxDistance_; }
		Float_t btag_jetBProbabilityBJetTags() const { return btag_jetBProbabilityBJetTags_; }
		Float_t btag_jetProbabilityBJetTags() const { return btag_jetProbabilityBJetTags_; }
		Float_t btag_trackCountingHighPurBJetTags() const { return btag_trackCountingHighPurBJetTags_; }
		Float_t btag_trackCountingHighEffBJetTags() const { return btag_trackCountingHighEffBJetTags_; }
		Float_t btag_simpleSecondaryVertexHighEffBJetTags() const { return btag_simpleSecondaryVertexHighEffBJetTags_; }
		Float_t btag_simpleSecondaryVertexHighPurBJetTags() const { return btag_simpleSecondaryVertexHighPurBJetTags_; }
		Float_t btag_combinedSecondaryVertexBJetTags() const { return btag_combinedSecondaryVertexBJetTags_; }
		Float_t btag_combinedSecondaryVertexRetrainedBJetTags() const { return btag_combinedSecondaryVertexRetrainedBJetTags_; }
		Float_t btag_combinedSecondaryVertexMVABJetTags() const { return btag_combinedSecondaryVertexMVABJetTags_; }
		Float_t btag_softMuonBJetTags() const { return btag_softMuonBJetTags_; }
		Float_t btag_softMuonByPtBJetTags() const { return btag_softMuonByPtBJetTags_; }
		Float_t	btag_softMuonByIP3dBJetTags() const { return btag_softMuonByIP3dBJetTags_; }
		Float_t	btag_softElectronByPtBJetTags() const { return btag_softElectronByPtBJetTags_; }
		Float_t	btag_softElectronByIP3dBJetTags() const { return btag_softElectronByIP3dBJetTags_; }
		Float_t	btag_combinedCSVJPBJetTags() const { return btag_combinedCSVJPBJetTags_; }
		Float_t	btag_combinedCSVJPSLBJetTags() const { return btag_combinedCSVJPSLBJetTags_; }
		Float_t	btag_combinedCSVSLBJetTags() const { return btag_combinedCSVSLBJetTags_; }
		Float_t	btag_softPFElectronRetrainedBJetsTags() const { return btag_softPFElectronRetrainedBJetsTags_; }
    Float_t	btag_softPFMuonRetrainedBJetsTags() const { return btag_softPFMuonRetrainedBJetsTags_; }

		std::map<std::string, float> getMistag_SF() const { 
			std::cout << mistag_SF_.size() << endl;
			return mistag_SF_;
		}
		std::map<std::string, float> getBtag_SF() const { return btag_SF_;}
		std::map<std::string, float> getMistag_SFerr() const { return mistag_SFerr_;}
		std::map<std::string, float> getBtag_SFerr() const { return btag_SFerr_;}


		Int_t partonFlavour() const {return partonFlavour_; }
		//Float_t partonFlavour() const {return partonFlavour_; }
		Bool_t isTopJet() const { return isTopJet_; }

		float getJetCorrFactor(std::string name) {

		  unsigned int size = sizeof(JetCorrName_)/sizeof(JetCorrName_[0]);

		  for (unsigned int i=0; i<size; i++) {

		    if ( JetCorrName_[i] == name )
		      return JetCorrValue_[i];

		  }

		  // if we reach this point, the correction factor was not found -> print all possible names
		  cout << "JetCorrFactor " << name << " was not found, possible names are: ";

		  for (unsigned int i=0; i<size; i++)
		    if (JetCorrName_[i] != "")
		      cout << JetCorrName_[i] << endl;

		  cout << endl;
		  
		  return 0;

		}

		virtual TString typeName() const { return "CatJet"; }

		void setJetType(Int_t jetType) { jetType_ = jetType; }
		void setNConstituents(Int_t nConstituents) { nConstituents_ = nConstituents; }
		void setJetArea(Float_t jetArea) { jetArea_ = jetArea; }
		void setMaxDistance(Float_t maxDistance) { maxDistance_ = maxDistance; }
		// btag
		void setBtag_jetBProbabilityBJetTags(Float_t btag_jetBProbabilityBJetTags) { btag_jetBProbabilityBJetTags_ = btag_jetBProbabilityBJetTags; }
		void setBtag_jetProbabilityBJetTags(Float_t btag_jetProbabilityBJetTags) { btag_jetProbabilityBJetTags_ = btag_jetProbabilityBJetTags; }
		void setBtag_trackCountingHighPurBJetTags(Float_t btag_trackCountingHighPurBJetTags) { btag_trackCountingHighPurBJetTags_ = btag_trackCountingHighPurBJetTags; }
		void setBtag_trackCountingHighEffBJetTags(Float_t btag_trackCountingHighEffBJetTags) { btag_trackCountingHighEffBJetTags_ = btag_trackCountingHighEffBJetTags; }
		void setBtag_simpleSecondaryVertexHighEffBJetTags(Float_t btag_simpleSecondaryVertexHighEffBJetTags) { btag_simpleSecondaryVertexHighEffBJetTags_ = btag_simpleSecondaryVertexHighEffBJetTags; }
		void setBtag_simpleSecondaryVertexHighPurBJetTags(Float_t btag_simpleSecondaryVertexHighPurBJetTags) { btag_simpleSecondaryVertexHighPurBJetTags_ = btag_simpleSecondaryVertexHighPurBJetTags; }
		void setBtag_combinedSecondaryVertexBJetTags(Float_t btag_combinedSecondaryVertexBJetTags) { btag_combinedSecondaryVertexBJetTags_ = btag_combinedSecondaryVertexBJetTags; }
		void setBtag_combinedSecondaryVertexRetrainedBJetTags(Float_t btag_combinedSecondaryVertexBJetTags) { btag_combinedSecondaryVertexRetrainedBJetTags_ = btag_combinedSecondaryVertexBJetTags; }
		void setBtag_combinedSecondaryVertexMVABJetTags(Float_t btag_combinedSecondaryVertexMVABJetTags) { btag_combinedSecondaryVertexMVABJetTags_ = btag_combinedSecondaryVertexMVABJetTags; }
		void setBtag_softMuonBJetTags(Float_t btag_softMuonBJetTags) { btag_softMuonBJetTags_ = btag_softMuonBJetTags; }
		void setBtag_softMuonByPtBJetTags(Float_t btag_softMuonByPtBJetTags) { btag_softMuonByPtBJetTags_ = btag_softMuonByPtBJetTags; }
		void setBtag_softMuonByIP3dBJetTags(Float_t btag_softMuonByIP3dBJetTags) { btag_softMuonByIP3dBJetTags_ = btag_softMuonByIP3dBJetTags; }
		void setBtag_softElectronByPtBJetTags(Float_t btag_softElectronByPtBJetTags) { btag_softElectronByPtBJetTags_ = btag_softElectronByPtBJetTags; }
		void setBtag_softElectronByIP3dBJetTags(Float_t btag_softElectronByIP3dBJetTags) { btag_softElectronByIP3dBJetTags_ = btag_softElectronByIP3dBJetTags; }
    void setBtag_combinedCSVJPBJetTags(Float_t btag_combinedCSVJPBJetTags) { btag_combinedCSVJPBJetTags_ = btag_combinedCSVJPBJetTags; }
    void setBtag_combinedCSVJPSLBJetTags(Float_t btag_combinedCSVJPSLBJetTags) { btag_combinedCSVJPSLBJetTags_ = btag_combinedCSVJPSLBJetTags; }
    void setBtag_combinedCSVSLBJetTags(Float_t btag_combinedCSVSLBJetTags) { btag_combinedCSVSLBJetTags_ = btag_combinedCSVSLBJetTags; }
    void setBtag_softPFElectronRetrainedBJetsTags(Float_t btag_softPFElectronRetrainedBJetsTags) { btag_softPFElectronRetrainedBJetsTags_ = btag_softPFElectronRetrainedBJetsTags; }
    void setBtag_softPFMuonRetrainedBJetsTags(Float_t btag_softPFMuonRetrainedBJetsTags) { btag_softPFMuonRetrainedBJetsTags_ = btag_softPFMuonRetrainedBJetsTags; }


		void setPartonFlavour(Int_t partonFlavour) { partonFlavour_ = partonFlavour; }
		void setIsTopJet(Bool_t isTopJet) { isTopJet_ = isTopJet; }

		//btag scalefactors
		void setMistag_SF(std::map<std::string, float> mistag_SF) { 
			for(std::map<std::string,float>::const_iterator it = mistag_SF.begin(); it != mistag_SF.end(); it++) {
				mistag_SF_[it->first] = it->second;
			}
		}
		void setBtag_SF(std::map<std::string, float> btag_SF) { 
			for(std::map<std::string,float>::const_iterator it = btag_SF.begin(); it != btag_SF.end(); it++) {
				btag_SF_[it->first] = it->second;
			}				
		}
		void setMistag_SFerr(std::map<std::string, float> mistag_SFerr) { 
			for(std::map<std::string,float>::const_iterator it = mistag_SFerr.begin(); it != mistag_SFerr.end(); it++) {
				mistag_SFerr_[it->first] = it->second;
			}				
		}
		void setBtag_SFerr(std::map<std::string, float> btag_SFerr) { 
			for(std::map<std::string,float>::const_iterator it = btag_SFerr.begin(); it != btag_SFerr.end(); it++) {
				btag_SFerr_[it->first] = it->second;
			}				
		}


		// JEC
		void setJetCorrFactor(int pos, std::string name, float factor)
		{  
		  JetCorrName_[pos] = name;
		  JetCorrValue_[pos] = factor;
		}

		friend std::ostream& operator<< (std::ostream& stream, const CatJet& jet)
		{
			stream << "CatJet - Charge=" << setw(2) << jet.charge() << " (Et,eta,phi)=("<< setw(10) << jet.Et() <<","<< setw(10) << jet.Eta() <<","<< setw(10) << jet.Phi() << ")"
				<< " vertex(x,y,z)=("<< jet.vx() <<","<< jet.vy() <<","<< jet.vz() << ")";
			return stream;
		};


	private:
		//Jet Info
		Int_t jetType_;                     // 0 = Unknown ; 1 = CaloJet ; 2 = PFJet
		Int_t nConstituents_;               // Number of constituents of the jet (calotowers for CaloJet / PFParticles for PFJet)
		Float_t jetArea_;                   // Jet area
		Float_t maxDistance_;               // Maximum distance from jet to constituent

		// jet correction factors
		std::string JetCorrName_[4]; 			// check in JetAnalyzer.cc that size is big enough to store all corrections!
		float JetCorrValue_[4];

		//btag Info
		Float_t btag_jetBProbabilityBJetTags_;
		Float_t btag_jetProbabilityBJetTags_;
		Float_t btag_trackCountingHighPurBJetTags_;
		Float_t btag_trackCountingHighEffBJetTags_;
		Float_t btag_simpleSecondaryVertexHighEffBJetTags_;
		Float_t btag_simpleSecondaryVertexHighPurBJetTags_;
		Float_t btag_combinedSecondaryVertexBJetTags_;
		Float_t btag_combinedSecondaryVertexRetrainedBJetTags_;
		Float_t btag_combinedSecondaryVertexMVABJetTags_;
		Float_t btag_softMuonBJetTags_;
		Float_t btag_softMuonByPtBJetTags_;
		Float_t btag_softMuonByIP3dBJetTags_;
		Float_t btag_softElectronByPtBJetTags_;
		Float_t btag_softElectronByIP3dBJetTags_;
		Float_t btag_combinedCSVJPBJetTags_;
		Float_t btag_combinedCSVJPSLBJetTags_;
		Float_t btag_combinedCSVSLBJetTags_;
		Float_t btag_softPFElectronRetrainedBJetsTags_;
		Float_t btag_softPFMuonRetrainedBJetsTags_;

		//btag scalefactors
		std::map<std::string, float> mistag_SF_;
		std::map<std::string, float> btag_SF_;
		std::map<std::string, float> mistag_SFerr_;
		std::map<std::string, float> btag_SFerr_;
		

		//MC info
		Int_t partonFlavour_;
		Bool_t isTopJet_;				// Is parton matched to the jet a decay product of the top quark ?

		ClassDef (CatJet,3);
	};
}

#endif
