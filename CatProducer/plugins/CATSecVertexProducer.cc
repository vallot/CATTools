#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CATTools/DataFormats/interface/SecVertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


using namespace edm;
using namespace std;
using namespace reco;


namespace cat {

  class CATSecVertexProducer : public edm::stream::EDProducer<> {
    public:
      explicit CATSecVertexProducer(const edm::ParameterSet & iConfig);

      void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

    private:
      void fitTransientTracks(SecVertexCollection* out_, reco::Vertex& goodPV, reco::TransientTrack& ttrack1, reco::TransientTrack& ttrack2, int pdgId);

      vector<cat::SecVertex> *out_;

      edm::EDGetTokenT<pat::PackedCandidateCollection>        pfSrc_;
      edm::EDGetTokenT<pat::MuonCollection>                   muonSrc_;
      edm::EDGetTokenT<pat::ElectronCollection>               elecSrc_;
      edm::EDGetTokenT<VertexCompositePtrCandidateCollection> secvtxSrc_;
      edm::EDGetTokenT<reco::TrackCollection>                 trackSrc_;
      edm::EDGetTokenT<reco::VertexCollection>                vertexLabel_;

      double rawMassMin_, rawMassMax_, massMin_, massMax_;

      double cut_minPt_, cut_maxEta_;
      double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
      int cut_trackNHit_, mode_, pfmode_;
      double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

      bool applyCuts_;

  };

} // namespace

cat::CATSecVertexProducer::CATSecVertexProducer(const edm::ParameterSet & iConfig) :
  pfSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  muonSrc_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  elecSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("elecSrc"))),
  secvtxSrc_(consumes<VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secvtxSrc"))),
  trackSrc_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
  produces<std::vector<cat::SecVertex> >();
  edm::ParameterSet trackPSet = iConfig.getParameter<edm::ParameterSet>("track");
  cut_minPt_ = trackPSet.getParameter<double>("minPt");
  cut_maxEta_ = trackPSet.getParameter<double>("maxEta");
  cut_trackChi2_ = trackPSet.getParameter<double>("chi2");
  cut_trackNHit_  = trackPSet.getParameter<int>("nHit");
  cut_trackSignif_ = trackPSet.getParameter<double>("signif");
  cut_DCA_ = trackPSet.getParameter<double>("DCA");

  edm::ParameterSet vertexPSet = iConfig.getParameter<edm::ParameterSet>("vertex");
  cut_vertexChi2_ = vertexPSet.getParameter<double>("chi2");
  cut_minLxy_ = vertexPSet.getParameter<double>("minLxy");
  cut_maxLxy_ = vertexPSet.getParameter<double>("maxLxy");
  cut_vtxSignif_ = vertexPSet.getParameter<double>("signif");

  rawMassMin_ = iConfig.getParameter<double>("rawMassMin");
  rawMassMax_ = iConfig.getParameter<double>("rawMassMax");
  massMin_ = iConfig.getParameter<double>("massMin");
  massMax_ = iConfig.getParameter<double>("massMax");

  applyCuts_ = iConfig.getParameter<bool>("applyCut");
  mode_ = iConfig.getParameter<int>("mode"); 
  pfmode_ = iConfig.getParameter<int>("pfmode"); 
}

  void
cat::CATSecVertexProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<pat::MuonCollection> muonSrc;
  Handle<pat::ElectronCollection> elecSrc;
  Handle<reco::TrackCollection> trackSrc;
  Handle<VertexCompositePtrCandidateCollection> secvtxSrc;
  Handle<pat::PackedCandidateCollection> pfSrc;
  if ( mode_ ==0 || mode_ == 1 ) {
    iEvent.getByToken(muonSrc_, muonSrc);
    iEvent.getByToken(elecSrc_, elecSrc);
  }
  if ( mode_== 1 ) {
    iEvent.getByToken(trackSrc_,trackSrc);
  }
  if ( mode_== 2 ) {
    iEvent.getByToken(secvtxSrc_, secvtxSrc);
  }
  if ( mode_ ==3 ) {
    iEvent.getByToken(pfSrc_, pfSrc);
  }

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);

  out_ = new std::vector<cat::SecVertex>();

  if ( recVtxs->empty() ) {
    unique_ptr<cat::SecVertexCollection> out(out_);
    iEvent.put(std::move(out));
    return;
  }
  reco::Vertex pv = recVtxs->at(0);

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);


  std::vector<TransientTrack> muTransTracks;
  std::vector<TransientTrack> elTransTracks;
  std::vector<TransientTrack> gTTracks;

  if ( mode_ == 0 || mode_ ==1 ) {
    for (const pat::Muon & aPatMuon : *muonSrc){
      reco::TrackRef trackRef;
      if ( aPatMuon.isGlobalMuon() ) trackRef = aPatMuon.globalTrack();
      else if ( aPatMuon.isTrackerMuon() ) trackRef = aPatMuon.innerTrack();
      else continue;
      reco::TransientTrack transTrack = trackBuilder->build(trackRef);
      muTransTracks.push_back( transTrack ) ;
    }
    for (const pat::Electron & aPatElectron : *elecSrc){
      if ( aPatElectron.gsfTrack().isNull() ) continue;
      reco::TransientTrack transTrack = trackBuilder->build(aPatElectron.gsfTrack());
      elTransTracks.push_back( transTrack ) ;
    }
  }

  if ( mode_==0 ) {
    if ( muTransTracks.size() >=2 ) { 
      for(auto ttrack1 = muTransTracks.begin() ; ttrack1 != muTransTracks.end()-1 ; ++ttrack1){
        for(auto ttrack2 = muTransTracks.begin()+(ttrack1-muTransTracks.begin())+1 ; ttrack2 != muTransTracks.end(); ++ttrack2) {
          fitTransientTracks(out_, pv, *ttrack1, *ttrack2, 13);
        }
      }
    }
    if( elTransTracks.size() >= 2) { 
      for(auto ttrack1 = elTransTracks.begin() ; ttrack1 != elTransTracks.end()-1 ; ++ttrack1){
        for(auto ttrack2 = elTransTracks.begin()+(ttrack1-elTransTracks.begin())+1 ; ttrack2 != elTransTracks.end(); ++ttrack2) {
          fitTransientTracks(out_, pv, *ttrack1, *ttrack2, 11);
        }
      }
    }
  }
  else if ( mode_==1 ) {
    for (auto it = trackSrc->begin(), ed = trackSrc->end(); it != ed; ++it) {
      unsigned int idx = it - trackSrc->begin();
      const auto  & aGeneralTrack = trackSrc->at(idx);

      reco::TransientTrack transTrack = trackBuilder->build(aGeneralTrack);
      if ( !transTrack.impactPointTSCP().isValid() ) continue;
      gTTracks.push_back(transTrack);
    }
    if ( muTransTracks.size() >0 ) {
      for(auto ttrack1 = muTransTracks.begin() ; ttrack1 != muTransTracks.end() ; ++ttrack1){
        for(auto ttrack2 = gTTracks.begin();     ttrack2 != gTTracks.end()      ; ++ttrack2) {
          fitTransientTracks(out_, pv, *ttrack1, *ttrack2, 13);
        }
      }
    }
    if ( elTransTracks.size() > 0 ) {
      for(auto ttrack1 = elTransTracks.begin() ; ttrack1 != elTransTracks.end() ; ++ttrack1){
        for(auto ttrack2 = gTTracks.begin() ; ttrack2 !=   gTTracks.end() ; ++ttrack2) {
          fitTransientTracks(out_, pv, *ttrack1, *ttrack2, 11);
        }
      }
    }
  }

  else if ( mode_ == 2) {
    for(VertexCompositePtrCandidateCollection::const_iterator sv = secvtxSrc->begin() ; sv != secvtxSrc->end() ; ++sv) {

      int pdgId = 0 ;
      // If the secvtx did not have 2 daughter(consist of 2 tracks), it must be ignored.
      if ( sv->numberOfDaughters() != 2 ) continue;
      const reco::Candidate *cand1, *cand2;
      cand1 = sv->daughter(0);
      cand2 = sv->daughter(1);
      if ( std::abs(cand1->pdgId()) == 211 && std::abs(cand2->pdgId()) == 211) continue;

      VertexCompositeCandidate aVC = VertexCompositeCandidate( sv->charge(), sv->p4(), sv->vertex(), sv->vertexCovariance(), sv->vertexChi2(), sv->vertexNdof(), pdgId, 0, true);
      cat::SecVertex aSecVtx(aVC);

      aSecVtx.setVProb(TMath::Prob( aVC.vertexChi2(), (int) aVC.vertexNdof()));
      auto vertex = aVC.vertex();
      GlobalPoint vtxPos( vertex.x(), vertex.y(), vertex.z() );

      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
      typedef ROOT::Math::SVector<double, 3> SVector3;

      SMatrixSym3D totalCov = pv.covariance() + aVC.vertexCovariance();
      SVector3 distanceVectorXY(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), 0.);

      double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
      double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;

      if(applyCuts_ && ( rVtxMag < cut_minLxy_ or rVtxMag > cut_maxLxy_ or rVtxMag / sigmaRvtxMag < cut_vtxSignif_ )) continue;

      SVector3 distanceVector3D(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), vertex.z() - pv.position().z());
      const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);

      aSecVtx.setLxy(rVtxMag);
      aSecVtx.setL3D(rVtxMag3D);
      //aSecVtx.setLeptonID(0,0);
      //aSecVtx.setTrackQuality(0,0);
      out_->push_back( aSecVtx); 
    }
  }
  else if ( mode_ == 3) {
    for(pat::PackedCandidateCollection::const_iterator cand1 = pfSrc->begin() ; cand1 != pfSrc->end()-1 ; ++cand1) {
      for(pat::PackedCandidateCollection::const_iterator cand2 = pfSrc->begin()+(cand1-pfSrc->begin()) +1 ; cand2 != pfSrc->end() ; ++cand2) {
        int pdgId = 0 ;
        if ( pfmode_ == 0 ) {
          if ( cand1->pdgId()*cand2->pdgId() == -121) pdgId = 11;
          if ( cand1->pdgId()*cand2->pdgId() == -169) pdgId = 13;
        }
        else if ( pfmode_ == 1 ) {
          if ( std::abs(cand1->pdgId()) == 11 || std::abs(cand2->pdgId()) == 11) pdgId = 11;
          if ( std::abs(cand1->pdgId()) == 13 || std::abs(cand2->pdgId()) == 13) pdgId = 13;
        }
        else if ( pfmode_ == 2) {
          //std::cout<<"cand1 mass : "<<cand1->mass()<<"  cand2 mass : "<<cand2->mass()<<std::endl;
          pdgId= 13; 
        }
        if (pdgId ==0 || cand1->charge()*cand2->charge() >= 0 ) continue; 

        auto track1 = cand1->pseudoTrack();
        auto track2 = cand2->pseudoTrack();
        auto transTrack1 = trackBuilder->build(track1);
        auto transTrack2 = trackBuilder->build(track2);
        try {
          fitTransientTracks(out_, pv, transTrack1, transTrack2, pdgId );
        } catch(std::exception& e) { std::cerr<<"Something error to fit TransientTracks : "<<e.what()<<std::endl; continue ; }
      }
    }
  }
  else { std::cerr<<"Can not found mode variable.Skip this event : "<<iEvent.id()<< std::endl; return ;}
  unique_ptr<cat::SecVertexCollection > out(out_);
  iEvent.put(std::move(out));
}

void cat::CATSecVertexProducer::fitTransientTracks(cat::SecVertexCollection* out_, reco::Vertex& goodPV, reco::TransientTrack& ttrack1, reco::TransientTrack& ttrack2, int pdgId)
{
  const double pvx = goodPV.position().x();
  const double pvy = goodPV.position().y();
  const double pvz = goodPV.position().z();
  double leptonMass = 0;
  if (pdgId == 11) leptonMass = 0.000511;
  if (pdgId == 13) leptonMass = 0.1056583715;

  FreeTrajectoryState ipStatePos ;
  FreeTrajectoryState ipStateNeg ;

  try {
    ipStatePos = ttrack1.impactPointTSCP().theState();
    ipStateNeg = ttrack2.impactPointTSCP().theState();
  } catch(std::exception& e) { std::cerr<<"theState()"<<std::endl; return ; }

  try{
    // Measure distance between tracks at their closest approach
    ClosestApproachInRPhi cApp;
    cApp.calculate(ipStatePos, ipStateNeg);
    if ( !cApp.status() ) return;
    const float dca = std::abs(cApp.distance());
    if ( applyCuts_ && (dca < 0. || dca > cut_DCA_ )) return;
    GlobalPoint cxPt = cApp.crossingPoint();
    if (applyCuts_ && (std::hypot(cxPt.x(), cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.)) return;
  } catch(std::exception& e) { std::cerr<<"closest approach"<<std::endl; return ; }

  // Build Vertex
  std::vector<TransientTrack> transTracks;
  transTracks.push_back(ttrack1 );
  transTracks.push_back(ttrack2 );
  KalmanVertexFitter fitter(true);
  TransientVertex transVertex;
  try{
    transVertex = fitter.vertex(transTracks);
  }catch(std::exception& e) { std::cerr<<"Kalman Vertex Fitting error for J/psi: "<<e.what()<<std::endl; return ; }

  if ( !transVertex.isValid() or transVertex.totalChiSquared() < 0. ) return;

  const reco::Vertex vertex = transVertex;
  if ( applyCuts_ && (vertex.normalizedChi2() > cut_vertexChi2_) ) return;

  std::vector<TransientTrack> refittedTracks;
  if ( transVertex.hasRefittedTracks() ) refittedTracks = transVertex.refittedTracks();

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z());
  SMatrixSym3D totalCov = goodPV.covariance() + vertex.covariance();
  SVector3 distanceVectorXY(vertex.x() - pvx, vertex.y() - pvy, 0.);

  double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
  double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;
  if(applyCuts_ && ( rVtxMag < cut_minLxy_ or rVtxMag > cut_maxLxy_ or rVtxMag / sigmaRvtxMag < cut_vtxSignif_ )) return;

  SVector3 distanceVector3D(vertex.x() - pvx, vertex.y() - pvy, vertex.z() - pvz);
  const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);

  // Cuts finished, now we create the candidates and push them back into the collections.
  std::unique_ptr<TrajectoryStateClosestToPoint> traj1;
  std::unique_ptr<TrajectoryStateClosestToPoint> traj2;

  try{
    if ( refittedTracks.empty() )
    {
      traj1.reset(new TrajectoryStateClosestToPoint(ttrack1.trajectoryStateClosestToPoint(vtxPos)));
      traj2.reset(new TrajectoryStateClosestToPoint(ttrack2.trajectoryStateClosestToPoint(vtxPos)));
    }
    else
    {
      TransientTrack* refTrack1 = 0, * refTrack2 = 0;
      for ( auto& refTrack : refittedTracks )
      {
        if ( refTrack.track().charge() < 0 ) refTrack1 = &refTrack;
        else refTrack2 = &refTrack;
      }
      if ( refTrack1 == 0 or refTrack2 == 0 ) return;
      traj1.reset(new TrajectoryStateClosestToPoint(refTrack1->trajectoryStateClosestToPoint(vtxPos)));
      traj2.reset(new TrajectoryStateClosestToPoint(refTrack2->trajectoryStateClosestToPoint(vtxPos)));
    }
    if( !traj1->isValid() or !traj2->isValid() ) return;
  } catch(std::exception& e) { std::cerr<<"TSCP errors : "<<e.what()<<std::endl; return ; }

  GlobalVector mom1(traj1->momentum());
  GlobalVector mom2(traj2->momentum());
  GlobalVector mom(mom1+mom2);

  //cleanup stuff we don't need anymore
  traj1.reset();
  traj2.reset();

  VertexCompositeCandidate* cand;
  const double candE1 = hypot(mom1.mag(), leptonMass);
  const double candE2 = hypot(mom2.mag(), leptonMass);
  double vtxChi2(vertex.chi2());
  double vtxNdof(vertex.ndof());
  try{
    Particle::Point vtx(vertex.x(), vertex.y(), vertex.z());
    const Vertex::CovarianceMatrix vtxCov(vertex.covariance());

    const math::XYZTLorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
    if ( applyCuts_ && (massMin_ > candLVec.mass() or massMax_ < candLVec.mass() )) return;

    // Match to muons
    cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
  }
  catch(std::exception& e) { std::cerr<<"Covariance Matrix error : "<<e.what()<<std::endl; return; }
  if (!cand) return;
  reco::LeafCandidate newLep1(+pdgId, math::XYZTLorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
  reco::LeafCandidate newLep2(-pdgId, math::XYZTLorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));

  cand->addDaughter(newLep1);
  cand->addDaughter(newLep2);

  cand->setPdgId(443);

  cat::SecVertex aSecVertex(*cand);
  aSecVertex.setVProb(TMath::Prob( vtxChi2, (int) vtxNdof));
  aSecVertex.setLxy(rVtxMag);
  aSecVertex.setL3D(rVtxMag3D);
  //aSecVertex.setLeptonID(0,0);
  //aSecVertex.setTrackQuality(0,0);

  out_->push_back(aSecVertex);

}
#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATSecVertexProducer);
