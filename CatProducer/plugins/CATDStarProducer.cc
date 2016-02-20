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

      const float gPionMass = 0.1396;
      const float gKaonMass = 0.4937;
      const float gD0Mass   = 1.86480;
      float d0MassWindow_;
      unsigned int maxNumPFCand_;
      bool applyCuts_;


  };
  bool pTComp( const reco::Candidate* a, const reco::Candidate* b) { return a->pt() > b->pt();  }
} // namespace

cat::CATDStarProducer::CATDStarProducer(const edm::ParameterSet & iConfig) :
  jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel")))
{
  produces<vector<cat::SecVertex> >("D0Cand");
  produces<vector<cat::SecVertex> >("DstarCand");

  maxNumPFCand_ = iConfig.getParameter<int>("maxNumPFCand");
  d0MassWindow_ = iConfig.getParameter<double>("d0MassWindow");
  applyCuts_ = iConfig.getParameter<bool>("applyCut");
}

  void
cat::CATDStarProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetSrc_, jetHandle);
  vector<cat::SecVertex>* D0_Out_    = new std::vector<cat::SecVertex>();
  vector<cat::SecVertex>* Dstar_Out_ = new std::vector<cat::SecVertex>();

  std::vector< std::vector<const reco::Candidate*> > jetsDaughters;

  for (const pat::Jet & aPatJet : *jetHandle){
    std::vector<const reco::Candidate*> jetDaughters;
    unsigned int dau_size = aPatJet.numberOfDaughters();
    if ( dau_size < 3 ) continue;
    if ( dau_size > maxNumPFCand_ ) break;
    for ( unsigned int pion_idx = 0 ; pion_idx< dau_size ; pion_idx++) {
      for ( unsigned int kaon_idx = 0 ; kaon_idx< dau_size ; kaon_idx++) {
        if ( pion_idx == kaon_idx ) continue;
        const reco::Candidate* pionCand = aPatJet.daughter(pion_idx);
        const reco::Candidate* kaonCand = aPatJet.daughter(kaon_idx);
        if ( abs(pionCand->pdgId()) != 211 || abs( kaonCand->pdgId()) != 211) continue;

        TLorentzVector pion, pion2, kaon, D0, Dstar ;
        pion.SetPtEtaPhiM(pionCand->pt(), pionCand->eta(), pionCand->phi(),gPionMass );
        kaon.SetPtEtaPhiM(kaonCand->pt(), kaonCand->eta(), kaonCand->phi(),gKaonMass );
        if ( pion.DeltaR(kaon) > 0.2 ) continue;

        D0 = pion+kaon;
        const math::XYZTLorentzVector lv( D0.Px(), D0.Py(), D0.Pz(), D0.E());
        VertexCompositeCandidate* D0Cand = new VertexCompositeCandidate(0, lv, Point(0,0,0), 521) ;  // + pdgId,
        D0Cand->addDaughter( *pionCand );
        D0Cand->addDaughter( *kaonCand );

        D0_Out_->push_back( cat::SecVertex(*D0Cand) );
        if ( abs( D0.M() - gD0Mass)  < d0MassWindow_ ) {
          for( unsigned int extra_pion_idx = 0 ;  extra_pion_idx < dau_size ; extra_pion_idx++) {
            if ( extra_pion_idx== pion_idx || extra_pion_idx == kaon_idx) continue;
            const reco::Candidate* pion2Cand = aPatJet.daughter(extra_pion_idx);
            if ( pion2Cand->pdgId() != 211) continue;
            pion2.SetPtEtaPhiM( pion2Cand->pt(), pion2Cand->eta(), pion2Cand->phi(), gPionMass);
            Dstar = pion+kaon+pion2;
            const math::XYZTLorentzVector lv2( Dstar.Px(), Dstar.Py(), Dstar.Pz(), Dstar.E());

            VertexCompositeCandidate* DstarCand = new VertexCompositeCandidate(pion2Cand->charge(), lv2, Point(0,0,0), 521211) ;  // + pdgId,
            DstarCand->addDaughter( *pionCand );
            DstarCand->addDaughter( *kaonCand );
            DstarCand->addDaughter( *pion2Cand );

            Dstar_Out_->push_back( cat::SecVertex(*DstarCand) );
          }
        }
      }

    }
  }
  auto_ptr<vector<cat::SecVertex> >    D0_Out(D0_Out_   );
  auto_ptr<vector<cat::SecVertex> > Dstar_Out(Dstar_Out_);
  iEvent.put(D0_Out   , "D0Cand");
  iEvent.put(Dstar_Out, "DstarCand");
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATDStarProducer);
