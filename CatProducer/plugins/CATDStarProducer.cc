#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <TLorentzVector.h>
#include "CATTools/DataFormats/interface/SecVertex.h"


#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include<memory>

using namespace edm;
using namespace std;
using namespace reco;

namespace cat {

  class CATDStarProducer : public edm::stream::EDProducer<> {
    public:
      explicit CATDStarProducer(const edm::ParameterSet & iConfig);
      virtual ~CATDStarProducer() { }

      void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
    private:

      edm::EDGetTokenT<edm::View<pat::Jet> >                    jetSrc_;
      edm::EDGetTokenT<reco::VertexCollection>                vertexLabel_;

      const float gPionMass = 0.1396;
      const float gKaonMass = 0.4937;
      const float gD0Mass   = 1.86480;
      float d0MassWindow_, maxDeltaR_ ,d0MassCut_;
      unsigned int maxNumPFCand_;
      bool applyCuts_;


  };
  //bool pTComp( const reco::Candidate* a, const reco::Candidate* b) { return a->pt() > b->pt();  }
} // namespace

cat::CATDStarProducer::CATDStarProducer(const edm::ParameterSet & iConfig) :
  jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel")))
{
  produces<vector<cat::SecVertex> >("D0Cand");
  produces<vector<cat::SecVertex> >("DstarCand");
  produces<vector<cat::SecVertex> >("JpsiCand");
  produces<vector<cat::SecVertex> >("JpsiMVACand");

  maxNumPFCand_ = iConfig.getParameter<int>("maxNumPFCand");
  d0MassWindow_ = iConfig.getParameter<double>("d0MassWindow");
  d0MassCut_ = iConfig.getParameter<double>("d0MassCut");
  maxDeltaR_  = iConfig.getParameter<double>("maxDeltaR");
  applyCuts_ = iConfig.getParameter<bool>("applyCut");
}
  void
cat::CATDStarProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  float dca;

  if ( recVtxs->empty() ) {
    auto_ptr<vector<cat::SecVertex> >    D0_Out(new vector<cat::SecVertex>());
    auto_ptr<vector<cat::SecVertex> > Dstar_Out(new std::vector<cat::SecVertex>());
    auto_ptr<vector<cat::SecVertex> >  Jpsi_Out(new std::vector<cat::SecVertex>());
    auto_ptr<vector<cat::SecVertex> >  JpsiMVA_Out(new std::vector<cat::SecVertex>());
    iEvent.put(D0_Out   , "D0Cand");
    iEvent.put(Dstar_Out, "DstarCand");
    iEvent.put(Jpsi_Out, "JpsiCand");
    iEvent.put(JpsiMVA_Out, "JpsiMVACand");
    return ; 
  }
  reco::Vertex pv = recVtxs->at(0);

  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetSrc_, jetHandle);

  auto_ptr<vector<cat::SecVertex> >    D0_Out(new vector<cat::SecVertex>());
  auto_ptr<vector<cat::SecVertex> > Dstar_Out(new std::vector<cat::SecVertex>());
  auto_ptr<vector<cat::SecVertex> >  Jpsi_Out(new std::vector<cat::SecVertex>());
  auto_ptr<vector<cat::SecVertex> >  JpsiMVA_Out(new std::vector<cat::SecVertex>());

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  typedef const pat::PackedCandidate ConstPC;
  //typedef std::shared_ptr<ConstPC> Shared_PCP;
  typedef ConstPC* Shared_PCP;

  for (const pat::Jet & aPatJet : *jetHandle){
    std::vector< Shared_PCP >  jetDaughters;
    std::vector<TransientTrack> tracks;
    unsigned int dau_size = aPatJet.numberOfDaughters();
    if ( dau_size < 3 ) continue;
    for( unsigned int idx = 0 ; idx < dau_size ; idx++) {
      jetDaughters.push_back( Shared_PCP(dynamic_cast<ConstPC*>(aPatJet.daughter(idx) ))); 

    }

    sort(jetDaughters.begin(), jetDaughters.end(), [](Shared_PCP a, Shared_PCP b) {return a->pt() > b->pt(); }); 

    if ( dau_size > maxNumPFCand_ ) dau_size = maxNumPFCand_;
    jetDaughters.resize( dau_size );


    for ( unsigned int lep1_idx = 0 ; lep1_idx< dau_size-1 ; lep1_idx++) {
      for ( unsigned int lep2_idx = lep1_idx+1 ; lep2_idx< dau_size ; lep2_idx++) {
        bool flag_jpsi = true;
        Shared_PCP lep1Cand = jetDaughters[lep1_idx];
        Shared_PCP lep2Cand = jetDaughters[lep2_idx];

        int pdgMul = lep1Cand->pdgId() * lep2Cand->pdgId();
        if ( pdgMul != -121 && pdgMul != -169 ) flag_jpsi = false;
        if ( !flag_jpsi && abs( pdgMul ) !=2321 && abs( pdgMul) != 2743 ) continue;


        //if ( reco::deltaR( *lep1Cand, *lep2Cand) > maxDeltaR_ ) continue;
        
        auto Jpsi = lep1Cand->p4()+ lep2Cand->p4();

        tracks.clear();
        KalmanVertexFitter fitter(true);
        TransientVertex t_vertex;
        reco::TransientTrack lep1Track, lep2Track;

        if ( lep1Cand->bestTrack() != nullptr && lep2Cand->bestTrack() != nullptr) {
          lep1Track = trackBuilder->build( lep1Cand->bestTrack());
          lep2Track = trackBuilder->build( lep2Cand->bestTrack());
          tracks.push_back( lep1Track);
          tracks.push_back( lep2Track);
        }
        else continue;
      
        try{
          t_vertex = fitter.vertex(tracks);
        }catch(std::exception& e) { std::cerr<<"Kalman Vertex Fitting error for Jpsi: "<<e.what()<<std::endl; }

        Point vx;
        reco::Vertex vertex;
        double vtxChi2=0.0;
        int vtxNdof=0;
        bool fit_jpsi = false; 
        if ( t_vertex.isValid() && t_vertex.totalChiSquared() > 0. )  {
          vertex = t_vertex;
          vx = Point(vertex.x(), vertex.y(), vertex.z());
          vtxChi2 = vertex.chi2(); 
          vtxNdof = (int)vertex.ndof();
          fit_jpsi = true;
        }
        else { vx = Point(0,0,0); continue; }


        const math::XYZTLorentzVector lv( Jpsi.px(), Jpsi.py(), Jpsi.pz(), Jpsi.E());
        auto vc = VertexCompositeCandidate(0, lv, vx, 443) ;  // + pdgId,
        cat::SecVertex JpsiCand(vc);
        JpsiCand.setJetDR( reco::deltaR( JpsiCand, aPatJet) );
        JpsiCand.setLegDR( reco::deltaR( *lep1Cand, *lep2Cand) );
        if ( fit_jpsi ) { 
          JpsiCand.setVProb( TMath::Prob( vtxChi2, vtxNdof));

          SVector3 distanceVectorXY(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), 0.);
          SVector3 distanceVector3D(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), vertex.z()- pv.position().z());
          double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
          double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
          JpsiCand.setLxy(rVtxMag);
          JpsiCand.setL3D(rVtxMag3D);
        }
        ClosestApproachInRPhi cApp;
        auto theLep1State = lep1Track.impactPointTSCP().theState();
        auto theLep2State = lep2Track.impactPointTSCP().theState();
        cApp.calculate(theLep1State, theLep2State);
        if ( cApp.status() ) { dca= std::abs(cApp.distance()); JpsiCand.set_dca(dca);}
        else { dca = -9 ; JpsiCand.set_dca(dca);}
        JpsiCand.set_dca(1,-9); 
        JpsiCand.set_dca(2,-9);
 
        JpsiCand.addDaughter( *lep1Cand );
        JpsiCand.addDaughter( *lep2Cand );

        JpsiCand.setTrackQuality( (int)lep1Cand->trackHighPurity(), (int)lep2Cand->trackHighPurity());
        if ( abs(lep1Cand->pdgId() ) == 13 || abs(lep2Cand->pdgId()) == 13 ) {
          int lep1ID = (int)lep1Cand->isStandAloneMuon() + (int)lep1Cand->isGlobalMuon()*2;
          int lep2ID = (int)lep2Cand->isStandAloneMuon() + (int)lep2Cand->isGlobalMuon()*2;
          JpsiCand.setLeptonID( lep1ID, lep2ID );
        }
        else JpsiCand.setLeptonID( -1, -1 );

        if( flag_jpsi ) Jpsi_Out->push_back( JpsiCand );
        JpsiMVA_Out->push_back( JpsiCand );

      }
    }

    for ( unsigned int pion_idx = 0 ; pion_idx< dau_size ; pion_idx++) {
      for ( unsigned int kaon_idx = 0 ; kaon_idx< dau_size ; kaon_idx++) {
        if ( pion_idx == kaon_idx ) continue;
        Shared_PCP pionCand = jetDaughters[pion_idx];
        std::shared_ptr<pat::PackedCandidate> kaonCand(jetDaughters[kaon_idx]->clone());
        kaonCand->setMass(gKaonMass);

        if ( abs(pionCand->pdgId()) != 211 || abs( kaonCand->pdgId()) != 211) continue;
        if ( pionCand->charge() * kaonCand->charge() != -1 ) continue;

        //if ( reco::deltaR( *pionCand, *kaonCand) > maxDeltaR_ ) continue;

        auto D0 = pionCand->p4()+ kaonCand->p4();
        if ( abs(D0.M() - gD0Mass) > d0MassCut_) continue;

        tracks.clear();
        TransientVertex t_vertex;
        reco::TransientTrack pionTrack, kaonTrack;
        KalmanVertexFitter fitter(true);
        if ( pionCand->bestTrack() != nullptr && kaonCand->bestTrack() != nullptr) {
          pionTrack = trackBuilder->build( pionCand->bestTrack());
          kaonTrack = trackBuilder->build( kaonCand->bestTrack());
          tracks.push_back( pionTrack);
          tracks.push_back( kaonTrack);
        }
        else continue;
      
        try{
          t_vertex = fitter.vertex(tracks);
        }catch(std::exception& e) { std::cerr<<"Kalman Vertex Fitting error for D0: "<<e.what()<<std::endl; }

        Point vx;
        reco::Vertex vertex;
        double vtxChi2=0.0;
        int vtxNdof=0;
        bool fit_d0 = false; 
        if ( t_vertex.isValid() && t_vertex.totalChiSquared() > 0. )  {
          vertex = t_vertex;
          vx = Point(vertex.x(), vertex.y(), vertex.z());
          vtxChi2 = vertex.chi2(); 
          vtxNdof = (int)vertex.ndof();
          fit_d0 = true;
          //printf(" D0 vertex => x : %e y: %e z: %e\n",vx.x(),vx.y(),vx.z());
        }
        else { vx = Point(0,0,0); continue; }


        const math::XYZTLorentzVector lv( D0.px(), D0.py(), D0.pz(), D0.E());
        auto vc = VertexCompositeCandidate(0, lv, vx, 421) ;  // + pdgId,
        cat::SecVertex D0Cand(vc);
        D0Cand.setJetDR( reco::deltaR( D0Cand, aPatJet) );
        D0Cand.setLegDR( reco::deltaR( *pionCand, *kaonCand) );
        if ( fit_d0) { 
          D0Cand.setVProb( TMath::Prob( vtxChi2, vtxNdof));

          SVector3 distanceVectorXY(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), 0.);
          SVector3 distanceVector3D(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), vertex.z()- pv.position().z());
          double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
          double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
          D0Cand.setLxy(rVtxMag);
          D0Cand.setL3D(rVtxMag3D);
        }
        ClosestApproachInRPhi cApp;
        auto thePionState = pionTrack.impactPointTSCP().theState();
        auto theKaonState = kaonTrack.impactPointTSCP().theState();
        cApp.calculate(thePionState, theKaonState);
        if ( cApp.status() ) { dca= std::abs(cApp.distance()); D0Cand.set_dca(dca);}
        else { dca = -9 ; D0Cand.set_dca(dca);}
        D0Cand.set_dca(1,-9); 
        D0Cand.set_dca(2,-9);
 
        D0Cand.addDaughter( *pionCand );
        D0Cand.addDaughter( *kaonCand );

        D0Cand.setTrackQuality( (int)pionCand->trackHighPurity(), (int)kaonCand->trackHighPurity());

        D0_Out->push_back( D0Cand );

        if ( abs( D0.M() - gD0Mass) < d0MassWindow_ ) {
          for( unsigned int extra_pion_idx = 0 ;  extra_pion_idx < dau_size ; extra_pion_idx++) {
            if ( extra_pion_idx== pion_idx || extra_pion_idx == kaon_idx) continue;
            Shared_PCP pion2Cand = jetDaughters[extra_pion_idx];
            if ( abs(pion2Cand->pdgId()) != 211) continue;
            //if ( reco::deltaR(D0Cand, *pion2Cand  )> maxDeltaR_) continue;
            auto Dstar = D0Cand.p4() + pion2Cand->p4();
            const math::XYZTLorentzVector lv2( Dstar.Px(), Dstar.Py(), Dstar.Pz(), Dstar.E());
            reco::TransientTrack pion2Track ; 
            if ( pion2Cand->bestTrack() != nullptr ) { 
              pion2Track = trackBuilder->build( pion2Cand->bestTrack());
            }
            else continue;
            tracks.clear();

            tracks.push_back( pionTrack);
            tracks.push_back( kaonTrack);
            tracks.push_back( pion2Track );
           
            bool fit_dstar = false; 
            try{
              t_vertex = fitter.vertex(tracks);
            }catch(std::exception& e) { std::cerr<<"Kalman Vertex Fitting error for D*: "<<e.what()<<std::endl; }
            if ( t_vertex.isValid() && t_vertex.totalChiSquared() > 0. )  {
              const reco::Vertex vertex = t_vertex; 
              vx = Point(vertex.x(), vertex.y(), vertex.z()); 
              vtxChi2 = vertex.chi2(); 
              vtxNdof = vertex.ndof();
              fit_dstar = true;
              //printf(" D* vertex => x : %e y: %e z: %e\n",vx.x(),vx.y(),vx.z());
            }
            else { vx = Point(0,0,0); continue; }
            
            auto vc2 = VertexCompositeCandidate(pion2Cand->charge(), lv2, vx, pion2Cand->charge()*413) ;  // + pdgId,
            cat::SecVertex DstarCand(vc2);
            DstarCand.setJetDR( reco::deltaR( DstarCand, aPatJet) );
            DstarCand.setLegDR( reco::deltaR(D0Cand, *pion2Cand ) );
            DstarCand.addDaughter( *pionCand );
            DstarCand.addDaughter( *kaonCand );
            DstarCand.addDaughter( *pion2Cand );
            DstarCand.setTrackQuality( (int)(pionCand->trackHighPurity()&kaonCand->trackHighPurity()), (int)pion2Cand->trackHighPurity());
            if ( fit_dstar) {
              DstarCand.setVProb( TMath::Prob( vtxChi2, (int) vtxNdof));

              SVector3 distanceVectorXY(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), 0.);
              SVector3 distanceVector3D(vertex.x() - pv.position().x(), vertex.y() - pv.position().y(), vertex.z()- pv.position().z());
              double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
              double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
              DstarCand.setLxy(rVtxMag);
              DstarCand.setL3D(rVtxMag3D);
            
            } 
            DstarCand.set_dca( 0, dca );
            auto thePion2State = pion2Track.impactPointTSCP().theState();
            cApp.calculate(thePionState, thePion2State);
            if ( cApp.status() ) DstarCand.set_dca(1, std::abs(cApp.distance()));
            else DstarCand.set_dca(1,-9);
            cApp.calculate(theKaonState, thePion2State);
            if ( cApp.status() ) DstarCand.set_dca(2, std::abs(cApp.distance()));
            else DstarCand.set_dca(2,-9);
      
            Dstar_Out->push_back( DstarCand );
          }
        }
      }

    }
  }
  iEvent.put(D0_Out   , "D0Cand");
  iEvent.put(Dstar_Out, "DstarCand");
  iEvent.put(Jpsi_Out, "JpsiCand");
  iEvent.put(JpsiMVA_Out, "JpsiMVACand");
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATDStarProducer);
