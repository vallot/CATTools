#ifndef CatMuon_h
#define CatMuon_h

#include "../interface/CatLepton.h"

using namespace std;

namespace cat
{
	class CatMuon : public CatLepton
	{

	public:

		CatMuon() :
			CatLepton()
			,vetoEm_(-9999.)
			,vetoHad_(-9999.)
			,chi2_(+9999.0)
		        ,nofValidHits_(-9999)
		        ,nofValidMuHits_(-9999)
		        ,nofValidPixelHits_(-9999)
			,nofMatchedStations_(-9999)
		        ,nofTrackerLayersWithMeasurement_(-9999)
			,algo_(-9999)
			,isPFMuon_(false)
			,id_(-9999)
			{;}

                CatMuon(const CatLepton& l) :
                        CatLepton(l)
                        ,vetoEm_(-9999.)
                        ,vetoHad_(-9999.)
                        ,chi2_(+9999.0)
                        ,nofValidHits_(-9999)
                        ,nofValidMuHits_(-9999)
                        ,nofValidPixelHits_(-9999)
                        ,nofMatchedStations_(-9999)
                        ,nofTrackerLayersWithMeasurement_(-9999)
                        ,algo_(-9999)
                        ,isPFMuon_(false)
                        ,id_(-9999)
                        {;}

		CatMuon(const CatMuon& muon) :
			CatLepton(muon)
			,vetoEm_(muon.vetoEm_)
			,vetoHad_(muon.vetoHad_)
			,chi2_(muon.chi2_)
			,nofValidHits_(muon.nofValidHits_)
		        ,nofValidMuHits_(muon.nofValidMuHits_)
		        ,nofValidPixelHits_(muon.nofValidPixelHits_)
			,nofMatchedStations_(muon.nofMatchedStations_)
		         ,nofTrackerLayersWithMeasurement_(muon.nofTrackerLayersWithMeasurement_)
			,algo_(muon.algo_)
			,isPFMuon_(muon.isPFMuon_)
			,id_(muon.id_)
			{;}

		CatMuon(Double_t px, Double_t py, Double_t pz, Double_t e) :
			CatLepton(px,py,pz,e)
			,vetoEm_(-9999.)
			,vetoHad_(-9999.)
			,chi2_(+9999.0)
			,nofValidHits_(-9999)
		        ,nofValidMuHits_(-9999)
		        ,nofValidPixelHits_(-9999)
			,nofMatchedStations_(-9999)
		        ,nofTrackerLayersWithMeasurement_(-9999)
			,algo_(-9999)
			,isPFMuon_(false)
			,id_(-9999)
			{;}

		CatMuon(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z) :
			CatLepton(px,py,pz,e,vtx_x,vtx_y,vtx_z)
			,vetoEm_(-9999.)
			,vetoHad_(-9999.)
			,chi2_(+9999.0)
			,nofValidHits_(-9999)
		        ,nofValidMuHits_(-9999)
		        ,nofValidPixelHits_(-9999)
			,nofMatchedStations_(-9999)
		        ,nofTrackerLayersWithMeasurement_(-9999)
			,algo_(-9999)
			,isPFMuon_(false)
			,id_(-9999)
			{;}

		CatMuon(Double_t px, Double_t py, Double_t pz, Double_t e, Double_t vtx_x, Double_t vtx_y, Double_t vtx_z, Int_t type, Int_t charge) :
			CatLepton(px,py,pz,e,vtx_x,vtx_y,vtx_z,type,charge)
			,vetoEm_(-9999.)
			,vetoHad_(-9999.)
			,chi2_(+9999.0)
			,nofValidHits_(-9999)
		        ,nofValidMuHits_(-9999)
		        ,nofValidPixelHits_(-9999)
			,nofMatchedStations_(-9999)
		        ,nofTrackerLayersWithMeasurement_(-9999)
			,algo_(-9999)
			,isPFMuon_(false)
			,id_(-9999)
			{;}

		CatMuon(const TLorentzVector &momentum) :
			CatLepton(momentum)
			,vetoEm_(-9999.)
			,vetoHad_(-9999.)
			,chi2_(+9999.0)
			,nofValidHits_(-9999)
		        ,nofValidMuHits_(-9999)
		        ,nofValidPixelHits_(-9999)
			,nofMatchedStations_(-9999)
		        ,nofTrackerLayersWithMeasurement_(-9999)
			,algo_(-9999)
			,isPFMuon_(false)
			,id_(-9999)
			{;}

		CatMuon(const TLorentzVector &momentum, const TVector3 &vertex, Int_t type, Int_t charge) :
			CatLepton(momentum, vertex, type, charge)
			,vetoEm_(-9999.)
			,vetoHad_(-9999.)
			,chi2_(+9999.0)
			,nofValidHits_(-9999)
		        ,nofValidMuHits_(-9999)
		        ,nofValidPixelHits_(-9999)
			,nofMatchedStations_(-9999)
		        ,nofTrackerLayersWithMeasurement_(-9999)
			,algo_(-9999)
			,isPFMuon_(false)
			,id_(-9999)
			{;}

		~CatMuon() {;}

		Int_t algo() const { return algo_;}
		Bool_t isGlobalMuon() const { return algo_ & 2; }
		Bool_t isTrackerMuon() const { return algo_ & 4; }
		Bool_t isStandAloneMuon() const { return algo_ & 8; }
		Bool_t isCaloMuon() const { return algo_ & 16; }
                Bool_t isPFMuon() const { return isPFMuon_; }
		Int_t id() const { return id_;}
		Bool_t idAllGlobalMuons() const { return id_ & 1; }
		Bool_t idAllTrackerMuons() const { return id_ & 2; }
		Bool_t idAllStandAloneMuons() const { return id_ & 4; }
		Bool_t idTrackerMuonArbitrated() const { return id_ & 8; }
		Bool_t idAllArbitrated() const { return id_ & 16; }
		Bool_t idGlobalMuonPromptTight() const { return id_ & 32; }
		Bool_t idTMLastStationLoose() const { return id_ & 64; }
		Bool_t idTMLastStationTight() const { return id_ & 128; }
		Bool_t idTMLastStationAngTight() const { return id_ & 256; }
		Bool_t idTMOneStationLoose() const { return id_ & 512; }
		Bool_t idTMOneStationTight() const { return id_ & 1024; }
		Bool_t idTMLastStationOptimizedLowPtLoose() const { return id_ & 2048; }
		Bool_t idTMLastStationOptimizedLowPtTight() const { return id_ & 4096; }
		Bool_t idTM2DCompatibilityLoose() const { return id_ & 8192; }
		Bool_t idTM2DCompatibilityTight() const { return id_ & 16384; }
		Float_t vetoEm() const { return vetoEm_;} 
		Float_t vetoHad() const { return vetoHad_;} 
		Float_t chi2() const { return chi2_;}
		Int_t nofValidHits() const { return nofValidHits_;} 
		Int_t nofValidMuHits() const { return nofValidMuHits_;}
		Int_t nofValidPixelHits() const { return nofValidPixelHits_;}
		Int_t nofMatchedStations() const { return nofMatchedStations_;}
		Int_t nofTrackerLayersWithMeasurement() const { return nofTrackerLayersWithMeasurement_;}
//		Float_t dB() const { return dB_; }
//		Float_t dBError() const { return dBError_; }
		virtual TString typeName() const { return "CatMuon"; }

		void setAlgo(Int_t algo) { algo_ = algo; }
                void setIsPFMuon(Bool_t isPFMuon) {isPFMuon_ = isPFMuon; }
		void setID(Int_t id) { id_ = id; }
		void setID(
			Int_t AllGlobalMuons,
			Int_t AllTrackerMuons,
			Int_t AllStandAloneMuons,
			Int_t TrackerMuonArbitrated,
			Int_t AllArbitrated,
			Int_t GlobalMuonPromptTight,
			Int_t TMLastStationLoose,
			Int_t TMLastStationTight,
			Int_t TMLastStationAngTight,
			Int_t TMOneStationLoose,
			Int_t TMOneStationTight,
			Int_t TMLastStationOptimizedLowPtLoose,
			Int_t TMLastStationOptimizedLowPtTight,
			Int_t TM2DCompatibilityLoose,
			Int_t TM2DCompatibilityTight)
		{
			id_ = AllGlobalMuons*1 + AllTrackerMuons*2 + AllStandAloneMuons*4 + TrackerMuonArbitrated*8 + AllArbitrated*16 + GlobalMuonPromptTight*32
				+ TMLastStationLoose*64 + TMLastStationTight*128 + TMLastStationAngTight*256 + TMOneStationLoose*512 + TMOneStationTight*1024 + TMLastStationOptimizedLowPtLoose*2048 
				+ TMLastStationOptimizedLowPtTight*4096 + TM2DCompatibilityLoose*8192 + TM2DCompatibilityTight*16384;
		}
		void setVetoEm(Float_t vetoEm) { vetoEm_ = vetoEm;}
		void setVetoHad(Float_t vetoHad) { vetoHad_ = vetoHad;}
		void setChi2(Float_t chi2){ chi2_ = chi2;}
		void setNofValidHits(Int_t nofValidHits){ nofValidHits_ = nofValidHits;}
		void setNofValidMuHits(Int_t x){ nofValidMuHits_ = x;}
		void setNofValidPixelHits(Int_t x){ nofValidPixelHits_ = x;}
		void setNofMatchedStations(Int_t x){ nofMatchedStations_ = x;}
		void setNofTrackerLayersWithMeasurement(Int_t x){ nofTrackerLayersWithMeasurement_ = x;}

		friend std::ostream& operator<< (std::ostream& stream, const CatMuon& muon)
		{
			stream << "CatMuon - Charge=" << muon.charge() << " (Et,eta,phi)=("<< muon.Et() <<","<< muon.Eta() <<","<< muon.Phi() << ")  vertex(x,y,z)=("<< muon.vx() <<","<< muon.vy() <<","<< muon.vz() << ")" << endl
				<< "            Type(G,T,S,C)=(" << muon.isGlobalMuon() << ","  << muon.isTrackerMuon() << ","  << muon.isStandAloneMuon() << ","  << muon.isCaloMuon() << ") "
				<< "  ID=(" << muon.idTrackerMuonArbitrated() << ","  << muon.idAllArbitrated() << ","  << muon.idGlobalMuonPromptTight() << ","  << muon.idTMLastStationLoose()
				<< ","  << muon.idTMLastStationTight() << ","  << muon.idTM2DCompatibilityLoose() << ","  << muon.idTM2DCompatibilityTight() << ")" << endl;
			return stream;
		};

	
	private:
		// Variables from reco::Muon
		Float_t vetoEm_;            // veto conesize is 0.07  in the ecal
		Float_t vetoHad_;           // veto conesize is 0.1  in the hcal
		Float_t chi2_;              // chi2 of global Muon
		Int_t nofValidHits_;        // nof hits of inner track
		Int_t nofValidMuHits_;      // nof hits on the global fit
		Int_t nofValidPixelHits_;   // nof pixel hits of inner track
		Int_t nofMatchedStations_;  // number of stations with matched segments
		Int_t nofTrackerLayersWithMeasurement_; 

		Int_t algo_; // binary => GlobalMuon=00010 , TrackerMuon=00100 , StandAloneMuon=01000 , CaloMuon=10000
		Bool_t isPFMuon_;
		Int_t id_; 		// MuonId coded in binary word id_ ==> TrackerMuonArbitrated=0000001 , AllArbitrated=0000010 , GlobalMuonPromptTight=0000100 ,
		// TMLastStationLoose=0001000 , TMLastStationTight=0010000 , TM2DCompatibilityLoose=0100000 , TM2DCompatibilityTight=1000000
		
		ClassDef (CatMuon,4);
	};
}

#endif
