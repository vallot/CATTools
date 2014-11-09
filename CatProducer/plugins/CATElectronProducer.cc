#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/isFinite.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATElectronProducer : public edm::EDProducer {
  public:
    explicit CATElectronProducer(const edm::ParameterSet & iConfig);
    virtual ~CATElectronProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    float getEffArea( float dR, float scEta );

    edm::EDGetTokenT<edm::View<pat::Electron> > src_;
    edm::EDGetTokenT<edm::View<pat::Vertex> > vertexLabel_;
    edm::EDGetTokenT<double> > rhoLabel_;
    bool runOnMC_;

  };

} // namespace

cat::CATElectronProducer::CATElectronProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src"))),
  vertexLabel_(consumes<edm::View<pat::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  rhoLabel_(consumes<double >(iConfig.getParameter<edm::InputTag>("rhoLabel"))),
  runOnMC_(iConfig.getParameter<bool>("runOnMC"))
{
  produces<std::vector<cat::Electron> >();
}

void 
cat::CATElectronProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
  Handle<View<pat::Electron> > src;
  iEvent.getByToken(src_, src);

  Handle<View<reco::Vertex> > recVtxs;
  iEvent.getByToken(vertexLabel_, recVtxs);
  reco::Vertex pv = recVtxs->at(0);

  Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rhoIso = std::max(*(rhoHandle.product()), 0.0);

  auto_ptr<vector<cat::Electron> >  out(new vector<cat::Electron>());

  for (View<pat::Electron>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
    unsigned int idx = it - src->begin();
    const pat::Electron & aPatElectron = src->at(idx);
    cat::Electron aElectron(aPatElectron);

    float scEta = aPatElectron.superCluster()->eta();
    double ecalpt = aPatElectron.ecalDrivenMomentum().pt();

    double elEffArea04 = getEffArea( 0.4, scEta);
    double chIso04 = aPatElectron.chargedHadronIso();
    double nhIso04 = aPatElectron.neutralHadronIso();
    double phIso04 = aPatElectron.photonIso();
    aElectron.setrelIso(0.4, chIso04, nhIso04, phIso04, elEffArea04, rhoIso, ecalpt);

    double elEffArea03 = getEffArea( 0.3, scEta);
    double chIso03 = aPatElectron.userIsolation("pat::User1Iso");
    double nhIso03 = aPatElectron.userIsolation("pat::User2Iso");
    double phIso03 = aPatElectron.userIsolation("pat::User3Iso");
    aElectron.setrelIso(0.3, chIso03, nhIso03, phIso03, elEffArea03, rhoIso, ecalpt);
    
    aElectron.setscEta( aPatElectron.superCluster()->eta());
    aElectron.setmvaTrigV0( aPatElectron.electronID("mvaTrigV0")) ;

    double dxy = fabs(aPatElectron.gsfTrack()->dxy(pv.position()));
    aElectron.setdxy( dxy ) ;
    double dz = fabs(aPatElectron.gsfTrack()->dz(pv.position()));
    aElectron.setdz( dz ) ;

    aElectron.setrho( rhoIso) ;
    
    aElectron.setconversionVeto( (aPatElectron.passConversionVeto() ) && (aPatElectron.gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=0) ) ;
    aElectron.setchargeIDFull( aPatElectron.isGsfCtfScPixChargeConsistent()) ;
    
    out->push_back(aElectron);


  }

  iEvent.put(out);
}

float 
cat::CATElectronProducer::getEffArea( float dR, float scEta) 
{
  ElectronEffectiveArea::ElectronEffectiveAreaTarget electronEATarget; 
  if ( runOnMC_ ) electronEATarget = ElectronEffectiveArea::kEleEAFall11MC;
  else electronEATarget = ElectronEffectiveArea::kEleEAData2012;

  if( dR < 0.35) 
    return ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, scEta, electronEATarget);
  else 
    return ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, scEta, electronEATarget);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATElectronProducer);
