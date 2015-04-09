#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/MuonReco/interface/MuonCocktails.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATMuonProducer : public edm::EDProducer {
  public:
    explicit CATMuonProducer(const edm::ParameterSet & iConfig);
    virtual ~CATMuonProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    bool mcMatch( const reco::Candidate::LorentzVector& lepton, Handle<reco::GenParticleCollection> genParticles );
    bool MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact );

  private:
    edm::EDGetTokenT<pat::MuonCollection> src_;
    edm::EDGetTokenT<pat::MuonCollection> shiftedEnDownSrc_;
    edm::EDGetTokenT<pat::MuonCollection> shiftedEnUpSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    edm::EDGetTokenT<reco::BeamSpot> beamLineSrc_;
    bool runOnMC_;
  };

} // namespace

cat::CATMuonProducer::CATMuonProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  shiftedEnDownSrc_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("shiftedEnDownSrc"))),
  shiftedEnUpSrc_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("shiftedEnUpSrc"))),
  mcLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  beamLineSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamLineSrc")))
{
  produces<std::vector<cat::Muon> >();
}

void 
cat::CATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
  runOnMC_ = !iEvent.isRealData();

  Handle<pat::MuonCollection> src;
  iEvent.getByToken(src_, src);

  Handle<reco::GenParticleCollection> genParticles;
  if (runOnMC_) iEvent.getByToken(mcLabel_,genParticles);
  
	// beamSpotHandle = beamSpot(in LQTreeMaker)  
  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamLineSrc_, beamSpotHandle);

	// recVtxs = primaryVertices(in LQTreeMaker)
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv= recVtxs->at(0);
   
  reco::BeamSpot beamSpot = *beamSpotHandle;
  reco::TrackBase::Point beamPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());

  edm::Handle<pat::MuonCollection> shiftedEnDownSrc;
  edm::Handle<pat::MuonCollection> shiftedEnUpSrc;
  if (runOnMC_){
    iEvent.getByToken(shiftedEnDownSrc_, shiftedEnDownSrc);
    iEvent.getByToken(shiftedEnUpSrc_, shiftedEnUpSrc);
  }

  auto_ptr<vector<cat::Muon> >  out(new vector<cat::Muon>());

	for (const pat::Muon & aPatMuon : *src) {
    cat::Muon aMuon(aPatMuon);

		aMuon.setEnergy( aPatMuon.energy() );   
 
		/// global muon ///
    if( aPatMuon.isGlobalMuon() ){
      aMuon.setGlobalCharge( aPatMuon.globalTrack()->charge() );
			aMuon.setGlobalEta( aPatMuon.globalTrack()->eta() );
			aMuon.setGlobalPt( aPatMuon.globalTrack()->pt() );
			aMuon.setGlobalPhi( aPatMuon.globalTrack()->phi() );
      aMuon.setPtError( aPatMuon.globalTrack()->ptError() );
      aMuon.setEtaError( aPatMuon.globalTrack()->etaError() );
			aMuon.setGlobalChi2( aPatMuon.globalTrack()->normalizedChi2() );
   		aMuon.setGlobalTrkValidHits( aPatMuon.globalTrack()->hitPattern().numberOfValidMuonHits() );
		}
    else{
			if( aPatMuon.track().isNonnull() && aPatMuon.track().isAvailable() ){
      	aMuon.setGlobalCharge( aPatMuon.track()->charge() );
				aMuon.setGlobalEta( aPatMuon.track()->eta() );
				aMuon.setGlobalPt( aPatMuon.track()->pt() );
				aMuon.setGlobalPhi( aPatMuon.track()->phi() );
      	aMuon.setPtError( aPatMuon.track()->ptError() );
      	aMuon.setEtaError( aPatMuon.track()->etaError() );
			}
			else{
				aMuon.setGlobalCharge( -999. );
				aMuon.setGlobalEta( -999. );
        aMuon.setPtError( -999. );
        aMuon.setEtaError( -999. );
			}
      aMuon.setGlobalChi2( -999. );
			aMuon.setGlobalTrkValidHits( -999. );
    }
   
    aMuon.setPFIsoR03ChargedHadron( aPatMuon.pfIsolationR03().sumChargedHadronPt );
    aMuon.setPFIsoR03NeutralHadron( aPatMuon.pfIsolationR03().sumNeutralHadronEt );
    aMuon.setPFIsoR03Photon( aPatMuon.pfIsolationR03().sumPhotonEt );
    aMuon.setPFIsoR04ChargedHadron( aPatMuon.pfIsolationR04().sumChargedHadronPt );
    aMuon.setPFIsoR04NeutralHadron( aPatMuon.pfIsolationR04().sumNeutralHadronEt );
    aMuon.setPFIsoR04Photon( aPatMuon.pfIsolationR04().sumPhotonEt );
    aMuon.setEcalVetoIso( aPatMuon.isolationR03().emVetoEt );
    aMuon.setHcalVetoIso( aPatMuon.isolationR03().hadVetoEt );
    aMuon.setPFIsoR03PU( aPatMuon.pfIsolationR03().sumPUPt );
    aMuon.setPFIsoR04PU( aPatMuon.pfIsolationR04().sumPUPt );
/*
		float genparPt = -999.;
		float genparEta= -999.;
		float genparPhi= -999.;
			
		if( runOnMC_ ){
			//aPatMuon.genParticleRefs().size() should be 0 or 1
			for(unsigned int igen = 0 ; igen < aPatMuon.genParticleRefs().size() ; ++igen ){
				genparPt = aPatMuon.genParticle(igen)->pt();
				genparEta= aPatMuon.genParticle(igen)->eta();
				genparPhi= aPatMuon.genParticle(igen)->phi();
			}
		}

		aMuon.setMatchedGenParticlePt( genparPt );
		aMuon.setMatchedGenParticleEta( genparEta );
		aMuon.setMatchedGenParticlePhi( genparPhi );
*/
		aMuon.setPrimaryVertexDXY( aPatMuon.dB() );
		aMuon.setPrimaryVertexDXYError( aPatMuon.edB() );

		float trkd0;
		if( aPatMuon.track().isNonnull() && aPatMuon.track().isAvailable() ){
			aMuon.setTrkVx( aPatMuon.track()->vx() );
    	aMuon.setTrkVy( aPatMuon.track()->vy() );
    	aMuon.setTrkVz( aPatMuon.track()->vz() );
			trkd0 = aPatMuon.track()->d0();
			if( beamSpotHandle.isValid() ){
				trkd0   = -(aPatMuon.track()->dxy( beamSpotHandle->position() ));
			}
			aMuon.setTrkD0Error( aPatMuon.track()->d0Error() );
			aMuon.setTrkPixelHits( aPatMuon.track()->hitPattern().numberOfValidPixelHits() );
			aMuon.setStationMatches( aPatMuon.numberOfMatchedStations() );
			aMuon.setTrackLayersWithMeasurement( aPatMuon.track()->hitPattern().trackerLayersWithMeasurement() );
			aMuon.setTrackerCharge( aPatMuon.track()->charge() );
		}
		else{
      aMuon.setTrkVx( -999. );
      aMuon.setTrkVy( -999. );
      aMuon.setTrkVz( -999. );
			trkd0 = -999. ;
			aMuon.setTrkD0Error( -999. );
		  aMuon.setTrkPixelHits( -999. );
      aMuon.setStationMatches( -999. );
      aMuon.setTrackLayersWithMeasurement( -999. );
			aMuon.setTrackerCharge( -999. );
		}
		aMuon.setTrkD0( trkd0 );
		aMuon.setIsPF( aPatMuon.isPFMuon() );
		aMuon.setIsGlobal( aPatMuon.isGlobalMuon() );
		aMuon.setIsTracker( aPatMuon.isTrackerMuon() );

		/// Cocktail Muon ///
		double bct_vtxDistXY_   = -9999., bct_vtxDistZ_    = -9999.;
		if( aPatMuon.isGlobalMuon() ){
			reco::TrackRef cocktail_track = (muon::tevOptimized(aPatMuon, 200, 17., 40., 0.25)).first; // this line gives error.
			aMuon.setCocktailPt( cocktail_track->pt() );
			aMuon.setCocktailEta( cocktail_track->eta() );
			aMuon.setCocktailPhi( cocktail_track->phi() );
			aMuon.setCocktailGlobalChi2( cocktail_track->normalizedChi2() );
			aMuon.setCocktailCharge( cocktail_track->charge() );
			if( recVtxs.isValid() ){
				double bct_bestdist3D = 999999.;
				for( reco::VertexCollection::const_iterator v_it=recVtxs->begin() ; v_it!=recVtxs->end() ; ++v_it ){
					double bct_distXY = cocktail_track->dxy(v_it->position());
					double bct_distZ  = cocktail_track->dz(v_it->position());
					double bct_dist3D = sqrt( pow(bct_distXY,2) + pow(bct_distZ,2) );
					if( bct_dist3D < bct_bestdist3D ){
					bct_bestdist3D = bct_dist3D;
					bct_vtxDistXY_   = bct_distXY;
					bct_vtxDistZ_    = bct_distZ;
					}
				}
			}
			aMuon.setCocktailTrkVtxDXY( bct_vtxDistXY_ );
			aMuon.setCocktailTrkVtxDZ( bct_vtxDistZ_ );
		}
		else{
      aMuon.setCocktailPt( -999. );
      aMuon.setCocktailEta( -999. );
      aMuon.setCocktailPhi( -999. );
      aMuon.setCocktailGlobalChi2( -999. );
      aMuon.setCocktailCharge( -999. );
      aMuon.setCocktailTrkVtxDXY( -999. );
      aMuon.setCocktailTrkVtxDZ( -999. );
		}

		if( aPatMuon.isStandAloneMuon() ){
			aMuon.setMuonSpecPt( aPatMuon.standAloneMuon()->pt() );
			aMuon.setMuonSpecEta( aPatMuon.standAloneMuon()->eta() );
			aMuon.setMuonSpecPhi( aPatMuon.standAloneMuon()->phi() );
			aMuon.setMuonSpecCharge( aPatMuon.standAloneMuon()->charge() );
		}
		else{
      aMuon.setMuonSpecPt( -999. );
      aMuon.setMuonSpecEta( -999. );
      aMuon.setMuonSpecPhi( -999. );
      aMuon.setMuonSpecCharge( -999. );
      aMuon.setMuonSpecE( -999. );
		}
		
		double minVtxDist3D = 9999.;
		double vtxDistXY_   = -9999.;
		double bt_minVtxDist3D = 9999.;
		int    bt_vtxIndex_    = -1;
		double bt_vtxDistXY_   = -9999.;
		double bt_vtxDistZ_    = -9999.;
		if( recVtxs.isValid() ){
			for( reco::VertexCollection::const_iterator v_it=recVtxs->begin() ; v_it!=recVtxs->end() ; ++v_it ){
				double distXY = -999., distZ = -999., dist3D = -999.;
				if( aPatMuon.track().isNonnull() && aPatMuon.track().isAvailable() ){
          distXY = aPatMuon.track()->dxy(v_it->position());
					distZ  = aPatMuon.track()->dz(v_it->position());
          dist3D = sqrt(pow(distXY,2) + pow(distZ,2));
				}
				if( dist3D < minVtxDist3D ){
        	minVtxDist3D = dist3D;
          vtxDistXY_   = distXY;
        }
        if( (aPatMuon.muonBestTrack()).isNonnull() ){
          double bt_distXY = aPatMuon.muonBestTrack()->dxy(v_it->position());
          double bt_distZ  = aPatMuon.muonBestTrack()->dz(v_it->position());
          double bt_dist3D = sqrt( pow(bt_distXY,2) + pow(bt_distZ,2) );
          if( bt_dist3D < bt_minVtxDist3D ){
            bt_minVtxDist3D = bt_dist3D;
            bt_vtxIndex_    = int(std::distance(recVtxs->begin(),v_it));
            bt_vtxDistXY_   = bt_distXY;
            bt_vtxDistZ_    = bt_distZ;
          }
        }	
			}
		}
		aMuon.setVtxDistXY( vtxDistXY_ );
		aMuon.setBestTrackVtxIndex( bt_vtxIndex_ );
		aMuon.setBestTrackVtxDistZ( bt_vtxDistZ_ );
		aMuon.setBestTrackVtxDistXY( bt_vtxDistXY_ );




		out->push_back(aMuon);

  }

  iEvent.put(out);
}

bool cat::CATMuonProducer::mcMatch( const reco::Candidate::LorentzVector& lepton, Handle<reco::GenParticleCollection> genParticles ){
  bool out = false;

  for (const reco::GenParticle & aGenPart : *genParticles){
    if( abs(aGenPart.pdgId()) != 13 ) continue;

    bool match = MatchObjects(lepton, aGenPart.p4(), false);

    if( match != true) continue;
   
    const reco::Candidate* mother = aGenPart.mother();
    while( mother != 0 ){
      if( abs(mother->pdgId()) == 23 || abs(mother->pdgId()) == 24 ) {
        out = true;
      }
      mother = mother->mother();
    }
  }

  return out;
}


bool cat::CATMuonProducer::MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact ) {
  double proEta = proObj.eta();
  double proPhi = proObj.phi();
  double proPt = proObj.pt();
  double pasEta = pasObj.eta();
  double pasPhi = pasObj.phi();
  double pasPt = pasObj.pt();

  double dRval = deltaR(proEta, proPhi, pasEta, pasPhi);
  double dPtRel = 999.0;
  if( proPt > 0.0 ) dPtRel = fabs( pasPt - proPt )/proPt;
  // If we are comparing two objects for which the candidates should
  // be exactly the same, cut hard. Otherwise take cuts from user.
  if( exact ) return ( dRval < 1e-3 && dPtRel < 1e-3 );
  else return ( dRval < 0.025 && dPtRel < 0.025 );
}



#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMuonProducer);
