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


using namespace edm;
using namespace std;

namespace cat {

  class CATMuonProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATMuonProducer(const edm::ParameterSet & iConfig);
    virtual ~CATMuonProducer() { }

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

    bool mcMatch( const reco::Candidate::LorentzVector& lepton, Handle<reco::GenParticleCollection> genParticles );
    bool MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact );
    double getMiniRelIso(edm::Handle<pat::PackedCandidateCollection> pfcands,  const reco::Candidate::LorentzVector& ptcl, double  r_iso_min, double r_iso_max , double kt_scale);
    
  private:
    edm::EDGetTokenT<pat::MuonCollection> src_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    edm::EDGetTokenT<pat::PackedCandidateCollection>        pfSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamLineSrc_;
    bool runOnMC_;

    typedef math::XYZPoint Point;
  };
} // namespace

cat::CATMuonProducer::CATMuonProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  mcLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  pfSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  beamLineSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamLineSrc")))
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





void
cat::CATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  Handle<pat::MuonCollection> src;
  iEvent.getByToken(src_, src);

  Handle<reco::GenParticleCollection> genParticles;
  if (runOnMC_) iEvent.getByToken(mcLabel_,genParticles);

  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamLineSrc_, beamSpotHandle);

  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);


  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv;
  if (recVtxs->size())
    pv = recVtxs->at(0);

  reco::BeamSpot beamSpot = *beamSpotHandle;
  reco::TrackBase::Point beamPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());

  GlobalPoint pVertex(pv.position().x(),pv.position().y(),pv.position().z());


  auto_ptr<vector<cat::Muon> >  out(new vector<cat::Muon>());
  for (const pat::Muon & aPatMuon : *src) {
    cat::Muon aMuon(aPatMuon);

    if (runOnMC_){
      aMuon.setGenParticleRef(aPatMuon.genParticleRef());
      aMuon.setMCMatched( mcMatch( aPatMuon.p4(), genParticles ) );      
    }

    aMuon.setChargedHadronIso04( aPatMuon.pfIsolationR04().sumChargedHadronPt );
    aMuon.setNeutralHadronIso04( aPatMuon.pfIsolationR04().sumNeutralHadronEt );
    aMuon.setPhotonIso04( aPatMuon.pfIsolationR04().sumPhotonEt );
    aMuon.setPUChargedHadronIso04( aPatMuon.pfIsolationR04().sumPUPt );

    aMuon.setChargedHadronIso03( aPatMuon.pfIsolationR03().sumChargedHadronPt );
    aMuon.setNeutralHadronIso03( aPatMuon.pfIsolationR03().sumNeutralHadronEt );
    aMuon.setPhotonIso03( aPatMuon.pfIsolationR03().sumPhotonEt );
    aMuon.setPUChargedHadronIso03( aPatMuon.pfIsolationR03().sumPUPt );



    //////////////// pfcands //////////////////         
    edm::Handle<pat::PackedCandidateCollection> pfcands;
    iEvent.getByToken(pfSrc_, pfcands);


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
  if( proPt > 0.0 ) dPtRel = std::abs( pasPt - proPt )/proPt;
  // If we are comparing two objects for which the candidates should
  // be exactly the same, cut hard. Otherwise take cuts from user.
  if( exact ) return ( dRval < 1e-3 && dPtRel < 1e-3 );
  else return ( dRval < 0.025 && dPtRel < 0.025 );
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATMuonProducer);
