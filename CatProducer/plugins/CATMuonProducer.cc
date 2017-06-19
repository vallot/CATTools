#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
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
#include "DataFormats/MuonReco/interface/MuonCocktails.h" // for cocktail muon
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "CATTools/CommonTools/interface/GenParticleHelper.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATMuonProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATMuonProducer(const edm::ParameterSet & iConfig);

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

    bool mcMatch( const reco::Candidate::LorentzVector& lepton, const std::vector<reco::GenParticleRef>& genParticles ) const;
    double getMiniRelIso(edm::Handle<pat::PackedCandidateCollection> pfcands,  const reco::Candidate::LorentzVector& ptcl, double  r_iso_min, double r_iso_max , double kt_scale);

  private:
    edm::EDGetTokenT<pat::MuonCollection> src_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    edm::EDGetTokenT<pat::PackedCandidateCollection>        pfSrc_;
    //edm::EDGetTokenT<reco::BeamSpot> beamLineSrc_;
    bool runOnMC_;

    typedef math::XYZPoint Point;
  };
} // namespace

cat::CATMuonProducer::CATMuonProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  mcLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  pfSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc")))
  //beamLineSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamLineSrc")))
{
  produces<std::vector<cat::Muon> >();
}



double cat::CATMuonProducer::getMiniRelIso(edm::Handle<pat::PackedCandidateCollection> pfcands,
					   const reco::Candidate::LorentzVector& ptcl,
					   double r_iso_min, double r_iso_max, double kt_scale){

  if (ptcl.pt()<5.) return 99999.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;

  double iso_nh(0.); double iso_ch(0.);
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh(0.5);
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl.pt()));
  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;
    double dr = deltaR(pfc, ptcl);
    if (dr > r_iso) continue;


    //////////////////  NEUTRALS  /////////////////////////
    if (pfc.charge()==0){
      if (pfc.pt()>ptThresh) {
        /////////// PHOTONS ////////////
        if (abs(pfc.pdgId())==22) {
          if(dr < deadcone_ph) continue;
          iso_ph += pfc.pt();
          /////////// NEUTRAL HADRONS ////////////
        } else if (abs(pfc.pdgId())==130) {
          if(dr < deadcone_nh) continue;
          iso_nh += pfc.pt();
        }
      }
      //////////////////  CHARGED from PV  /////////////////////////
    } else if (pfc.fromPV()>1){
      if (abs(pfc.pdgId())==211) {
        if(dr < deadcone_ch) continue;
        iso_ch += pfc.pt();
      }
      //////////////////  CHARGED from PU  /////////////////////////
    } else {
      if (pfc.pt()>ptThresh){
        if(dr < deadcone_pu) continue;
        iso_pu += pfc.pt();
      }
    }
  }
  double iso(0.);
  iso = iso_ph + iso_nh;
  iso -= 0.5*iso_pu;
  if (iso>0) iso += iso_ch;
  else iso = iso_ch;

  iso = iso/ptcl.pt();

  return iso;
}


void cat::CATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  edm::Handle<pat::MuonCollection> src;
  iEvent.getByToken(src_, src);

  edm::Handle<reco::GenParticleCollection> genParticles;
  std::vector<reco::GenParticleRef> genMuons;
  if (runOnMC_) {
    iEvent.getByToken(mcLabel_,genParticles);
    genMuons = cat::selectGenParticles(genParticles, 13, {23,24});
  }

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  /*
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamLineSrc_, beamSpotHandle);
  const reco::BeamSpot beamSpot = *beamSpotHandle;
  const reco::TrackBase::Point beamPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());
  */

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv;
  if ( !recVtxs->empty() ) pv = recVtxs->at(0);
  const GlobalPoint pVertex(pv.position().x(),pv.position().y(),pv.position().z());

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfSrc_, pfcands);

  std::unique_ptr<cat::MuonCollection>  out(new cat::MuonCollection());
  for (const pat::Muon & aPatMuon : *src) {
    cat::Muon aMuon(aPatMuon);

    if (runOnMC_){
      aMuon.setGenParticleRef(aPatMuon.genParticleRef());
      aMuon.setMCMatched( mcMatch(aPatMuon.p4(), genMuons) );
    }

    aMuon.setChargedHadronIso04( aPatMuon.pfIsolationR04().sumChargedHadronPt );
    aMuon.setNeutralHadronIso04( aPatMuon.pfIsolationR04().sumNeutralHadronEt );
    aMuon.setPhotonIso04( aPatMuon.pfIsolationR04().sumPhotonEt );
    aMuon.setPUChargedHadronIso04( aPatMuon.pfIsolationR04().sumPUPt );

    aMuon.setChargedHadronIso03( aPatMuon.pfIsolationR03().sumChargedHadronPt );
    aMuon.setNeutralHadronIso03( aPatMuon.pfIsolationR03().sumNeutralHadronEt );
    aMuon.setPhotonIso03( aPatMuon.pfIsolationR03().sumPhotonEt );
    aMuon.setPUChargedHadronIso03( aPatMuon.pfIsolationR03().sumPUPt );

    aMuon.setMiniRelIso(getMiniRelIso( pfcands, aMuon.p4(), 0.05, 0.2, 10.));
    aMuon.setIsGlobalMuon( aPatMuon.isGlobalMuon() );
    aMuon.setIsTrackerMuon( aPatMuon.isTrackerMuon() );
    aMuon.setIsPF( aPatMuon.isPFMuon() );
    aMuon.setIsTight( aPatMuon.isTightMuon(pv) );
    aMuon.setIsMedium( aPatMuon.isMediumMuon() );
    aMuon.setIsLoose( aPatMuon.isLooseMuon() );
    aMuon.setIsSoftMuon( aPatMuon.isSoftMuon(pv) );

    aMuon.setNumberOfMatchedStations( aPatMuon.numberOfMatchedStations() );
    aMuon.setNumberOfValidHits( aPatMuon.numberOfValidHits() );

    if ( aPatMuon.globalTrack().isNonnull() && aPatMuon.globalTrack().isAvailable() ){
      aMuon.setNormalizedChi2( aPatMuon.normChi2() );
      aMuon.setNumberOfValidMuonHits( aPatMuon.globalTrack()->hitPattern().numberOfValidMuonHits() );
    }

    if ( aPatMuon.innerTrack().isNonnull() && aPatMuon.innerTrack().isAvailable() ){
      aMuon.setNumberOfValidPixelHits( aPatMuon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      aMuon.setTackerLayersWithMeasurement( aPatMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement() );
    }

    aMuon.setDxy( aPatMuon.muonBestTrack()->dxy(pv.position()) );
    aMuon.setDz( aPatMuon.muonBestTrack()->dz(pv.position()) );
    aMuon.setVertex(Point(aPatMuon.muonBestTrack()->vx(),aPatMuon.muonBestTrack()->vy(),aPatMuon.muonBestTrack()->vz()));


    reco::TrackRef  mutrack = aPatMuon.get<reco::TrackRef>();
    if (mutrack.isNull()){
      mutrack=aPatMuon.get<reco::TrackRef,reco::StandAloneMuonTag>();
    }
    reco::TransientTrack mutranstrack = trackBuilder->build( mutrack ) ;

    TrajectoryStateOnSurface  muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());

    if (muTSOS.isValid()){
      std::pair<bool,Measurement1D> muIPpair = IPTools::signedTransverseImpactParameter(mutranstrack, muTSOS.globalDirection(),pv);
      float muSignificanceIP = muIPpair.second.significance();
      aMuon.setIpSignficance(muSignificanceIP);
    }

    out->push_back(aMuon);
  }

  iEvent.put(std::move(out));
}

bool cat::CATMuonProducer::mcMatch( const reco::Candidate::LorentzVector& lepton, const std::vector<reco::GenParticleRef>& genMuons ) const
{
  for (const auto& aGenPart : genMuons){
    if ( cat::isMatchedByDRDPt(lepton, aGenPart->p4(), false) ) return true;
  }
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMuonProducer);
