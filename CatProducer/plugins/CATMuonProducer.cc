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
#include "DataFormats/MuonReco/interface/MuonCocktails.h" // for cocktail muon

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
    
  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamLineSrc_, beamSpotHandle);

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
  int j = 0;
  for (const pat::Muon & aPatMuon : *src) {
    cat::Muon aMuon(aPatMuon);

    /// cat default variables ///
    
    if (runOnMC_){
      aMuon.setShiftedEnDown(shiftedEnDownSrc->at(j).pt() );
      aMuon.setShiftedEnUp(shiftedEnUpSrc->at(j).pt() );
    }
    ++j;

    double chIso04 = aPatMuon.pfIsolationR04().sumChargedHadronPt;
    double nhIso04 = aPatMuon.pfIsolationR04().sumNeutralHadronEt;
    double phIso04 = aPatMuon.pfIsolationR04().sumPhotonEt;
    double puIso04 = aPatMuon.pfIsolationR04().sumPUPt;
    aMuon.setChargedHadronIso04( chIso04 );
    aMuon.setNeutralHadronIso04( nhIso04 );
    aMuon.setPhotonIso04( phIso04 );
    aMuon.setPUChargedHadronIso04( puIso04 );

    double chIso03 = aPatMuon.pfIsolationR03().sumChargedHadronPt;
    double nhIso03 = aPatMuon.pfIsolationR03().sumNeutralHadronEt;
    double phIso03 = aPatMuon.pfIsolationR03().sumPhotonEt;
    double puIso03 = aPatMuon.pfIsolationR03().sumPUPt;
    aMuon.setChargedHadronIso03( chIso03 );
    aMuon.setNeutralHadronIso03( nhIso03 );
    aMuon.setPhotonIso03( phIso03 );
    aMuon.setPUChargedHadronIso03( puIso03 );
    
    aMuon.setIsGlobalMuon( aPatMuon.isGlobalMuon() );
    aMuon.setIsPFMuon( aPatMuon.isPFMuon() );
    aMuon.setIsSoftMuon( aPatMuon.isSoftMuon(pv) );

    if (runOnMC_){
      aMuon.setGenParticleRef(aPatMuon.genParticleRef());
      aMuon.setMCMatched( mcMatch( aPatMuon.p4(), genParticles ) );
    }
    
    aMuon.setNumberOfMatchedStations( aPatMuon.numberOfMatchedStations() );

    if ( aPatMuon.globalTrack().isNonnull() && aPatMuon.globalTrack().isAvailable() ) {
      aMuon.setNormalizedChi2( aPatMuon.normChi2() );
      aMuon.setNumberOfValidMuonHits( aPatMuon.globalTrack()->hitPattern().numberOfValidMuonHits() );
    }

    if ( aPatMuon.innerTrack().isNonnull() && aPatMuon.innerTrack().isAvailable() ){
      aMuon.setNumberOfValidHits( aPatMuon.numberOfValidHits() );
      aMuon.setNumberOfValidPixelHits( aPatMuon.innerTrack()->hitPattern().numberOfValidPixelHits() );
      aMuon.setTackerLayersWithMeasurement( aPatMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement() ); 
    }
    double dxy = aPatMuon.muonBestTrack()->dxy(pv.position()); // fabs() removed
    aMuon.setDxy( dxy );
    double dz = aPatMuon.muonBestTrack()->dz(pv.position()); // fabs() removed
    aMuon.setDz( dz ); 
 
    /// SKTree variables ///
    
    aMuon.setEcalVetoIso( aPatMuon.isolationR03().emVetoEt );
    aMuon.setHcalVetoIso( aPatMuon.isolationR03().hadVetoEt );
    if( aPatMuon.track().isNonnull() && aPatMuon.track().isAvailable() ){
      aMuon.setTrkVx( aPatMuon.track()->vx() );
      aMuon.setTrkVy( aPatMuon.track()->vy() );
      aMuon.setTrkVz( aPatMuon.track()->vz() );
    }
    else{
      aMuon.setTrkVx( -999. );
      aMuon.setTrkVy( -999. );
      aMuon.setTrkVz( -999. );
    }
    aMuon.setIsTracker( aPatMuon.isTrackerMuon() );
    
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

