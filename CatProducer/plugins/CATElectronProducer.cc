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

using namespace edm;
using namespace std;

namespace cat {

  class CATElectronProducer : public edm::EDProducer {
  public:
    explicit CATElectronProducer(const edm::ParameterSet & iConfig);
    virtual ~CATElectronProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:
    edm::InputTag src_;
    edm::InputTag vertexLabel_;
    edm::InputTag rhoLabel_;

  };

} // namespace

cat::CATElectronProducer::CATElectronProducer(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>( "src" )),
  vertexLabel_(iConfig.getParameter<edm::InputTag>( "vertexLabel" )),
  rhoLabel_(iConfig.getParameter<edm::InputTag>("rhoLabel"))
{
  produces<std::vector<cat::Electron> >();
}

void 
cat::CATElectronProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{
  Handle<View<pat::Electron> > src;
  iEvent.getByLabel(src_, src);

  Handle<View<reco::Vertex> > recVtxs;
  iEvent.getByLabel(vertexLabel_, recVtxs);
  reco::Vertex pv = recVtxs->at(0);

  Handle<double> rhoHandle;
  iEvent.getByLabel(rhoLabel_, rhoHandle);
  const double rho = *(rhoHandle.product());



  auto_ptr<vector<cat::Electron> >  out(new vector<cat::Electron>());

  for (View<pat::Electron>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
    unsigned int idx = it - src->begin();
    const pat::Electron & aPatElectron = src->at(idx);
    cat::Electron aElectron(aPatElectron);

    aElectron.setChargedHadronIso04( aPatElectron.chargedHadronIso() );
    aElectron.setNeutralHadronIso04( aPatElectron.neutralHadronIso() );
    aElectron.setPhotonIso04( aPatElectron.photonIso() );
    aElectron.setPUChargedHadronIso04( aPatElectron.puChargedHadronIso() );

    aElectron.setChargedHadronIso03( aPatElectron.userIsolation("pat::User1Iso") );
    aElectron.setNeutralHadronIso03( aPatElectron.userIsolation("pat::User2Iso") );
    aElectron.setPhotonIso03( aPatElectron.userIsolation("pat::User3Iso") );
    aElectron.setPUChargedHadronIso03( aPatElectron.userIsolation("pat::User4Iso") );
    
    aElectron.setscEta( aPatElectron.superCluster()->eta());
    aElectron.setmva( aPatElectron.electronID("mvaTrigV0")) ;

    double dxy = fabs(aPatElectron.gsfTrack()->dxy(pv.position()));
    aElectron.setdxy( dxy ) ;
    double dz = fabs(aPatElectron.gsfTrack()->dz(pv.position()));
    aElectron.setdz( dz ) ;

    aElectron.setrho( rho) ;

    aElectron.setconversionVeto( (aPatElectron.passConversionVeto() ) && (aPatElectron.gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=0) ) ;
    aElectron.setchargeIDFull( aPatElectron.isGsfCtfScPixChargeConsistent()) ;
    

    out->push_back(aElectron);
  }

  iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATElectronProducer);
