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
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATMuonProducer : public edm::EDProducer {
  public:
    explicit CATMuonProducer(const edm::ParameterSet & iConfig);
    virtual ~CATMuonProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    bool mcMatch( const reco::Candidate::LorentzVector& lepton, const edm::Handle<edm::View<reco::GenParticle> > & genParticles );
    bool MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact );

  private:
    edm::InputTag src_;
    edm::InputTag mcLabel_;
    edm::InputTag vertexLabel_;
    edm::InputTag beamLineSrc_;

  };

} // namespace

cat::CATMuonProducer::CATMuonProducer(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>( "src" )),
  mcLabel_(iConfig.getParameter<edm::InputTag>( "mcLabel" )),
  vertexLabel_(iConfig.getParameter<edm::InputTag>( "vertexLabel" )),
  beamLineSrc_(iConfig.getParameter<edm::InputTag>( "beamLineSrc" ))
{
  produces<std::vector<cat::Muon> >();
}

void 
cat::CATMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
  Handle<View<pat::Muon> > src;
  iEvent.getByLabel(src_, src);
 
  Handle<View<reco::GenParticle> > genParticles;
  iEvent.getByLabel(mcLabel_,genParticles);
    
  Handle<View<reco::Vertex> > recVtxs;
  iEvent.getByLabel(vertexLabel_,recVtxs);

  Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(beamLineSrc_, beamSpotHandle);

  reco::Vertex pv = recVtxs->at(0);
   
  reco::BeamSpot beamSpot = *beamSpotHandle;
  reco::TrackBase::Point beamPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());
  // //  beamPoint = reco::TrackBase::Point ( beamSpot.x0(), beamSpot.y0(), beamSpot.z0() );  
 
  auto_ptr<vector<cat::Muon> >  out(new vector<cat::Muon>());

  for (View<pat::Muon>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
    unsigned int idx = it - src->begin();
    const pat::Muon & aPatMuon = src->at(idx);

    //    bool mcMatched = mcMatch( aPatMuon.p4(), genParticles );

    cat::Muon aMuon(aPatMuon);

    aMuon.setChargedHadronIso04( aPatMuon.chargedHadronIso() );
    aMuon.setNeutralHadronIso04( aPatMuon.neutralHadronIso() );
    aMuon.setPhotonIso04( aPatMuon.photonIso() );
    aMuon.setPUChargedHadronIso04( aPatMuon.puChargedHadronIso() );

    // aMuon.setChargedHadronIso03( aPatMuon.userIsolation("pat::User1Iso") );
    // aMuon.setNeutralHadronIso03( aPatMuon.userIsolation("pat::User2Iso") );
    // aMuon.setPhotonIso03( aPatMuon.userIsolation("pat::User3Iso") );
    // aMuon.setPUChargedHadronIso03( aPatMuon.userIsolation("pat::User4Iso") );

    aMuon.setIsGlobalMuon( aPatMuon.isGlobalMuon() );
    aMuon.setIsPFMuon( aPatMuon.isPFMuon() );
    aMuon.setIsTightMuon( aPatMuon.isTightMuon(pv) );
    aMuon.setIsLooseMuon( aPatMuon.isLooseMuon() );
    aMuon.setIsSoftMuon( aPatMuon.isSoftMuon(pv) );

    // aMuon.setMCMatched( mcMatched );
    
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
    double dxy = fabs(aPatMuon.muonBestTrack()->dxy(pv.position()));
    aMuon.setDxy( dxy );
    double dz = fabs(aPatMuon.muonBestTrack()->dz(pv.position()));
    aMuon.setDz( dz ); 
 
    out->push_back(aMuon);
  }

  iEvent.put(out);
}

bool cat::CATMuonProducer::mcMatch( const reco::Candidate::LorentzVector& lepton, const edm::Handle<edm::View<reco::GenParticle> > & genParticles ){

  bool out = false;

  for (edm::View<reco::GenParticle>::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {
    int genId = mcIter->pdgId();

    if( abs(genId) != 13 ) continue;

    bool match = MatchObjects(lepton, mcIter->p4(), false);

    if( match != true) continue;
   
    const reco::Candidate* mother = mcIter->mother();
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
