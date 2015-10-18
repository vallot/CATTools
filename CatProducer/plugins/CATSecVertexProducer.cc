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
    virtual ~CATSecVertexProducer() { }

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    void fitTransientTracks(reco::Vertex goodPV, std::vector<TransientTrack> transTracks, int pdgId) const;

    vector<cat::SecVertex> *out_;

    edm::EDGetTokenT<pat::MuonCollection> muonSrc_;
    edm::EDGetTokenT<pat::ElectronCollection> elecSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;

    double rawMassMin_, rawMassMax_, massMin_, massMax_;

    double cut_minPt_, cut_maxEta_;
    double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
    int cut_trackNHit_;
    double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;
  };

} // namespace

cat::CATSecVertexProducer::CATSecVertexProducer(const edm::ParameterSet & iConfig) :
  muonSrc_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  elecSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("elecSrc"))),
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
}

void
cat::CATSecVertexProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<pat::MuonCollection> muonSrc;
  iEvent.getByToken(muonSrc_, muonSrc);

  Handle<pat::ElectronCollection> elecSrc;
  iEvent.getByToken(elecSrc_, elecSrc);

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);

  out_ = new std::vector<cat::SecVertex>();

  if ( recVtxs->empty() || (muonSrc->size() < 2 && elecSrc->size() < 2))
    return;

  reco::Vertex pv = recVtxs->at(0);

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  //make muon transient tracks
  std::vector<TransientTrack> muTransTracks;
  for (const pat::Muon & aPatMuon : *muonSrc){
    reco::TrackRef trackRef;
    if ( aPatMuon.isGlobalMuon() ) trackRef = aPatMuon.globalTrack();
    else if ( aPatMuon.isTrackerMuon() ) trackRef = aPatMuon.innerTrack();
    else continue;

    reco::TransientTrack transTrack = trackBuilder->build(trackRef);
    if ( !transTrack.impactPointTSCP().isValid() ) continue;

    muTransTracks.push_back(transTrack);
  }
  fitTransientTracks(pv, muTransTracks, 13);

  //make electron transient tracks
  std::vector<TransientTrack> elTransTracks;
  for (const pat::Electron & aPatElectron : *elecSrc){
    if ( aPatElectron.gsfTrack().isNull() ) continue;

    reco::TransientTrack transTrack = trackBuilder->build(aPatElectron.gsfTrack());
    if ( !transTrack.impactPointTSCP().isValid() ) continue;

    elTransTracks.push_back(transTrack);
  }
  fitTransientTracks(pv, elTransTracks, 11);

  auto_ptr<vector<cat::SecVertex> > out(out_);
  iEvent.put(out);
}

void
cat::CATSecVertexProducer::fitTransientTracks(reco::Vertex goodPV, std::vector<TransientTrack> transTracks, int pdgId) const
{
  const double pvx = goodPV.position().x();
  const double pvy = goodPV.position().y();
  const double pvz = goodPV.position().z();
  double leptonMass = 0;
  if (pdgId == 11) leptonMass = 0.000511;
  if (pdgId == 13) leptonMass = 0.1056583715;
  int ipos=-1,ineg=-1;

  for ( auto& transTrackPos : transTracks )
    {
      ++ipos;
      if (transTrackPos.charge() < 0) continue;
      FreeTrajectoryState ipStatePos = transTrackPos.impactPointTSCP().theState();
      ineg=-1;
      for ( auto& transTrackNeg : transTracks )
	{
	  ++ineg;
	  if (transTrackNeg.charge() > 0) continue;
	  FreeTrajectoryState ipStateNeg = transTrackNeg.impactPointTSCP().theState();

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(ipStatePos, ipStateNeg);
	  if ( !cApp.status() ) continue;
	  const float dca = std::abs(cApp.distance());
	  if ( dca < 0. || dca > cut_DCA_ ) continue;
	  GlobalPoint cxPt = cApp.crossingPoint();
	  if (std::hypot(cxPt.x(), cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue;
	  //TrajectoryStateClosestToPoint caState1 = transTrackPos.trajectoryStateClosestToPoint(cxPt);
	  //TrajectoryStateClosestToPoint caState2 = transTrackNeg.trajectoryStateClosestToPoint(cxPt);
	  //if ( !caState1.isValid() or !caState2.isValid() ) continue;

	  // Build Vertex
	  std::vector<TransientTrack> transTracks;
	  transTracks.push_back(transTrackPos);
	  transTracks.push_back(transTrackNeg);
	  KalmanVertexFitter fitter(true);
	  TransientVertex transVertex = fitter.vertex(transTracks);

	  if ( !transVertex.isValid() or transVertex.totalChiSquared() < 0. ) continue;

	  const reco::Vertex vertex = transVertex;
	  if ( vertex.normalizedChi2() > cut_vertexChi2_ ) continue;

	  std::vector<TransientTrack> refittedTracks;
	  if ( transVertex.hasRefittedTracks() ) refittedTracks = transVertex.refittedTracks();

	  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
	  typedef ROOT::Math::SVector<double, 3> SVector3;

	  GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z());
	  SMatrixSym3D totalCov = goodPV.covariance() + vertex.covariance();
	  SVector3 distanceVectorXY(vertex.x() - pvx, vertex.y() - pvy, 0.);

	  double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
	  double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;
	  if( rVtxMag < cut_minLxy_ or rVtxMag > cut_maxLxy_ or rVtxMag / sigmaRvtxMag < cut_vtxSignif_ ) continue;

	  SVector3 distanceVector3D(vertex.x() - pvx, vertex.y() - pvy, vertex.z() - pvz);
	  const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);

	  // Cuts finished, now we create the candidates and push them back into the collections.
	  std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
	  std::auto_ptr<TrajectoryStateClosestToPoint> traj2;

	  if ( refittedTracks.empty() )
	    {
	      traj1.reset(new TrajectoryStateClosestToPoint(transTrackPos.trajectoryStateClosestToPoint(vtxPos)));
	      traj2.reset(new TrajectoryStateClosestToPoint(transTrackNeg.trajectoryStateClosestToPoint(vtxPos)));
	    }
	  else
	    {
	      TransientTrack* refTrack1 = 0, * refTrack2 = 0;
	      for ( auto& refTrack : refittedTracks )
		{
		  if ( refTrack.track().charge() < 0 ) refTrack1 = &refTrack;
		  else refTrack2 = &refTrack;
		}
	      if ( refTrack1 == 0 or refTrack2 == 0 ) continue;
	      traj1.reset(new TrajectoryStateClosestToPoint(refTrack1->trajectoryStateClosestToPoint(vtxPos)));
	      traj2.reset(new TrajectoryStateClosestToPoint(refTrack2->trajectoryStateClosestToPoint(vtxPos)));
	    }
	  if( !traj1->isValid() or !traj2->isValid() ) continue;

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
	  if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;

	  // Match to muons
	  VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
	  if (!cand) continue;
	  reco::LeafCandidate newLep1(+pdgId, math::XYZTLorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
	  reco::LeafCandidate newLep2(-pdgId, math::XYZTLorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));

	  cand->addDaughter(newLep1);
	  cand->addDaughter(newLep2);

	  cand->setPdgId(pdgId);

	  cat::SecVertex aSecVertex(*cand);
	  aSecVertex.setVProb(TMath::Prob( vtxChi2, (int) vtxNdof));
	  aSecVertex.setLxy(rVtxMag);
	  aSecVertex.setL3D(rVtxMag3D);
	  aSecVertex.setInts(ipos,ineg);

	  out_->push_back(aSecVertex);
	}
    }

}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATSecVertexProducer);
