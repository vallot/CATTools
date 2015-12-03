#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CATTools/DataFormats/interface/SecVertexV2.h"
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

  class CATSecVertexV4Producer : public edm::stream::EDProducer<> {
    public:
      explicit CATSecVertexV4Producer(const edm::ParameterSet & iConfig);
      virtual ~CATSecVertexV4Producer() { }

      void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

    private:
      void fitTransientTracks(SecVertexV2Collection* out_, reco::Vertex& goodPV, reco::TransientTrack& ttrack1, reco::TransientTrack& ttrack2, int pdgId);

      vector<cat::SecVertexV2> *out_;

      //edm::EDGetTokenT<pat::MuonCollection>     muonSrc_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfSrc_;
      edm::EDGetTokenT<reco::VertexCollection>  vertexLabel_;

      double rawMassMin_, rawMassMax_, massMin_, massMax_;

      double cut_minPt_, cut_maxEta_;
      double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
      int cut_trackNHit_;
      double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

      bool applyCuts_;

  };

} // namespace

cat::CATSecVertexV4Producer::CATSecVertexV4Producer(const edm::ParameterSet & iConfig) :
  pfSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
  produces<std::vector<cat::SecVertexV2> >();
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
}

  void
cat::CATSecVertexV4Producer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<pat::PackedCandidateCollection> pfSrc;
  iEvent.getByToken(pfSrc_, pfSrc);

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);

  out_ = new std::vector<cat::SecVertexV2>();

  if ( recVtxs->empty() ) {
    auto_ptr<cat::SecVertexV2Collection> out(out_);
    iEvent.put(out);
    return;
  }
  reco::Vertex pv = recVtxs->at(0);

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  for(pat::PackedCandidateCollection::const_iterator cand1 = pfSrc->begin() ; cand1 != pfSrc->end()-1 ; ++cand1) {
    for(pat::PackedCandidateCollection::const_iterator cand2 = pfSrc->begin()+(cand1-pfSrc->begin()) +1 ; cand2 != pfSrc->end() ; ++cand2) {
    int pdgId = 0 ;
    if ( std::abs(cand1->pdgId()) == 11 || std::abs(cand2->pdgId()) == 11) pdgId = 11;
    if ( std::abs(cand1->pdgId()) == 13 || std::abs(cand2->pdgId()) == 13) pdgId = 13;

    if (pdgId ==0 || cand1->charge()*cand2->charge() >= 0 ) continue; 

    try {
      auto track1 = cand1->pseudoTrack();
      auto track2 = cand2->pseudoTrack();
      auto transTrack1 = trackBuilder->build(track1);
      auto transTrack2 = trackBuilder->build(track2);
      fitTransientTracks(out_, pv, transTrack1, transTrack2, pdgId );
    } catch(std::exception& e) { std::cerr<<"ll"<<std::endl; continue ; }
    }
  }
  auto_ptr<cat::SecVertexV2Collection > out(out_);
  iEvent.put(out);
}

void cat::CATSecVertexV4Producer::fitTransientTracks(cat::SecVertexV2Collection* out_, reco::Vertex& goodPV, reco::TransientTrack& ttrack1, reco::TransientTrack& ttrack2, int pdgId)
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
  //TrajectoryStateClosestToPoint caState1 = leptonTRack.trajectoryStateClosestToPoint(cxPt);
  //TrajectoryStateClosestToPoint caState2 = transTrackNeg.trajectoryStateClosestToPoint(cxPt);
  //if ( !caState1.isValid() or !caState2.isValid() ) continue;
  } catch(std::exception& e) { std::cerr<<"closest approach"<<std::endl; return ; }

  // Build Vertex
  std::vector<TransientTrack> transTracks;
  transTracks.push_back(ttrack1 );
  transTracks.push_back(ttrack2 );
  KalmanVertexFitter fitter(true);
  TransientVertex transVertex = fitter.vertex(transTracks);

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
  std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
  std::auto_ptr<TrajectoryStateClosestToPoint> traj2;

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

  GlobalVector mom1(traj1->momentum());
  GlobalVector mom2(traj2->momentum());
  GlobalVector mom(mom1+mom2);

  //cleanup stuff we don't need anymore
  traj1.reset();
  traj2.reset();

  Particle::Point vtx(vertex.x(), vertex.y(), vertex.z());
  const Vertex::CovarianceMatrix vtxCov(vertex.covariance());
  double vtxChi2(vertex.chi2());
  double vtxNdof(vertex.ndof());

  const double candE1 = hypot(mom1.mag(), leptonMass);
  const double candE2 = hypot(mom2.mag(), leptonMass);

  const math::XYZTLorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
  if ( applyCuts_ && (massMin_ > candLVec.mass() or massMax_ < candLVec.mass() )) return;

  // Match to muons
  VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
  if (!cand) return;
  reco::LeafCandidate newLep1(+pdgId, math::XYZTLorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
  reco::LeafCandidate newLep2(-pdgId, math::XYZTLorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));

  cand->addDaughter(newLep1);
  cand->addDaughter(newLep2);

  cand->setPdgId(pdgId);

  cat::SecVertexV2 aSecVertex(*cand);
  aSecVertex.setVProb(TMath::Prob( vtxChi2, (int) vtxNdof));
  aSecVertex.setLxy(rVtxMag);
  aSecVertex.setL3D(rVtxMag3D);
  aSecVertex.setInts(0,0);

  out_->push_back(aSecVertex);

}
#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATSecVertexV4Producer);
