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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "FWCore/Utilities/interface/transform.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h" //

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATElectronProducer : public edm::EDProducer {
  public:
    explicit CATElectronProducer(const edm::ParameterSet & iConfig);
    virtual ~CATElectronProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
    bool mcMatch( const reco::Candidate::LorentzVector& lepton, const edm::Handle<reco::GenParticleCollection> & genParticles );
    bool MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact );

  private:
    float getEffArea( float dR, float scEta );

    edm::EDGetTokenT<edm::View<pat::Electron> > src_;
    edm::EDGetTokenT<edm::View<pat::Electron> > shiftedEnDownSrc_;
    edm::EDGetTokenT<edm::View<pat::Electron> > shiftedEnUpSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
    edm::EDGetTokenT<reco::BeamSpot> beamLineSrc_;
    edm::EDGetTokenT<double> rhoLabel_;
    bool runOnMC_;
 
    typedef std::pair<std::string, edm::InputTag> NameTag;
    std::vector<NameTag> elecIDSrcs_;
    std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > elecIDTokens_;
  };

} // namespace

cat::CATElectronProducer::CATElectronProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src"))),
  shiftedEnDownSrc_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("shiftedEnDownSrc"))),
  shiftedEnUpSrc_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("shiftedEnUpSrc"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  mcLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  beamLineSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamLineSrc"))),
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel")))
{
  produces<std::vector<cat::Electron> >();
  if (iConfig.existsAs<edm::ParameterSet>("electronIDSources")) {
    edm::ParameterSet idps = iConfig.getParameter<edm::ParameterSet>("electronIDSources");
    std::vector<std::string> names = idps.getParameterNamesForType<edm::InputTag>();
    for (std::vector<std::string>::const_iterator it = names.begin(), ed = names.end(); it != ed; ++it) {
      auto inputTag = idps.getParameter<edm::InputTag>(*it);
      elecIDSrcs_.push_back(NameTag(inputTag.instance(), inputTag));
    }
    elecIDTokens_ = edm::vector_transform(elecIDSrcs_, [this](NameTag const & tag){return mayConsume<edm::ValueMap<bool> >(tag.second);});
  }
}

void 
cat::CATElectronProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;

  runOnMC_ = !iEvent.isRealData();

  Handle<edm::View<pat::Electron> > src;
  iEvent.getByToken(src_, src);

  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(mcLabel_,genParticles);

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_, recVtxs);
  reco::Vertex pv= recVtxs->at(0);

  Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rhoIso = std::max(*(rhoHandle.product()), 0.0);

  edm::Handle<edm::View<pat::Electron> > shiftedEnDownSrc;
  edm::Handle<edm::View<pat::Electron> > shiftedEnUpSrc;
  if (runOnMC_){
    iEvent.getByToken(shiftedEnDownSrc_, shiftedEnDownSrc);
    iEvent.getByToken(shiftedEnUpSrc_, shiftedEnUpSrc);
  }
  
  std::vector<edm::Handle<edm::ValueMap<bool> > > idhandles;
  std::vector<pat::Electron::IdPair>               ids;
  idhandles.resize(elecIDSrcs_.size());
  ids.resize(elecIDSrcs_.size());
  for (size_t i = 0; i < elecIDSrcs_.size(); ++i) {
    iEvent.getByToken(elecIDTokens_[i], idhandles[i]);
    ids[i].first = elecIDSrcs_[i].first;
  }

  auto_ptr<vector<cat::Electron> >  out(new vector<cat::Electron>());
  int j = 0;
  for (const pat::Electron &aPatElectron : *src){
    cat::Electron aElectron(aPatElectron);
    auto elecsRef = src->refAt(j);

    if (runOnMC_){
      aElectron.setShiftedEnDown(shiftedEnDownSrc->at(j).pt() );
      aElectron.setShiftedEnUp(shiftedEnUpSrc->at(j).pt() );
      aElectron.setGenParticleRef(aPatElectron.genParticleRef());
    }

    bool mcMatched = mcMatch( aPatElectron.p4(), genParticles );
    aElectron.setMCMatched( mcMatched );

    aElectron.setChargedHadronIso04( aPatElectron.chargedHadronIso() );
    aElectron.setNeutralHadronIso04( aPatElectron.neutralHadronIso() );
    aElectron.setPhotonIso04( aPatElectron.photonIso() );
    aElectron.setPUChargedHadronIso04( aPatElectron.puChargedHadronIso() );

    aElectron.setChargedHadronIso03( aPatElectron.pfIsolationVariables().sumChargedHadronPt );
    aElectron.setNeutralHadronIso03( aPatElectron.pfIsolationVariables().sumNeutralHadronEt );
    aElectron.setPhotonIso03( aPatElectron.pfIsolationVariables().sumPhotonEt );
    aElectron.setPUChargedHadronIso03( aPatElectron.pfIsolationVariables().sumPUPt );
 
    float scEta = aPatElectron.superCluster()->eta();
    double ecalpt = aPatElectron.ecalDrivenMomentum().pt();

    double elEffArea04 = getEffArea( 0.4, scEta);
    double chIso04 = aElectron.chargedHadronIso(0.4);
    double nhIso04 = aElectron.neutralHadronIso(0.4);
    double phIso04 = aElectron.photonIso(0.4);
    aElectron.setrelIso(0.4, chIso04, nhIso04, phIso04, elEffArea04, rhoIso, ecalpt);

    double elEffArea03 = getEffArea( 0.3, scEta);
    double chIso03 = aElectron.chargedHadronIso(0.3);
    double nhIso03 = aElectron.neutralHadronIso(0.3);
    double phIso03 = aElectron.photonIso(0.3);
    aElectron.setrelIso(0.3, chIso03, nhIso03, phIso03, elEffArea03, rhoIso, ecalpt);

    aElectron.setIsEB( aPatElectron.isEB());
    aElectron.setIsEE( aPatElectron.isEE());
    aElectron.setTrackerDrivenSeed( aPatElectron.trackerDrivenSeed());
    aElectron.setEcalDrivenSeed( aPatElectron.ecalDrivenSeed());

    aElectron.setTrkIsoDR03( aPatElectron.dr03TkSumPt());
    aElectron.setEcalIsoDR03( aPatElectron.dr03EcalRecHitSumEt());
    aElectron.setHcalIsoDR03( aPatElectron.dr03HcalTowerSumEt());
    aElectron.setTrkIso( aPatElectron.trackIso());
    aElectron.setEcalIso( aPatElectron.ecalIso());
    aElectron.setHcalIso( aPatElectron.hcalIso());

    aElectron.setDeltaPhiTrkSC( aPatElectron.deltaPhiSuperClusterTrackAtVtx());
    aElectron.setDeltaEtaTrkSC( aPatElectron.deltaEtaSuperClusterTrackAtVtx());

    aElectron.setSigmaIEtaIEta( aPatElectron.sigmaIetaIeta());

    aElectron.setHoE( aPatElectron.hadronicOverEm());
    aElectron.setCaloEnergy( aPatElectron.caloEnergy());
    aElectron.setESuperClusterOverP( aPatElectron.eSuperClusterOverP());
    aElectron.setTrackVx( aPatElectron.gsfTrack()->vx());
    aElectron.setTrackVy( aPatElectron.gsfTrack()->vy());
    aElectron.setTrackVz( aPatElectron.gsfTrack()->vz());

    aElectron.setNumberOfBrems( aPatElectron.numberOfBrems());
    aElectron.setFbrem( aPatElectron.fbrem());
    aElectron.setPrimaryVertexDXY( fabs(aPatElectron.dB()));
    aElectron.setPrimaryVertexDXYError( fabs(aPatElectron.edB()));
    aElectron.setTrackPt( aPatElectron.gsfTrack()->pt());
    aElectron.setTrackValidFractionOfHits( aPatElectron.gsfTrack()->validFraction());    
 
    float minVtxDist3D = 9999.;
    int vertexIndex_ = -1;
    float vertexDistXY_ = -9999.;
    float vertexDistZ_ = -9999.;

    aElectron.setPrimaryVertexDXY( fabs(aPatElectron.dB()));
    aElectron.setPrimaryVertexDXYError( fabs(aPatElectron.edB()));
    aElectron.setTrackPt( aPatElectron.gsfTrack()->pt());
    aElectron.setTrackValidFractionOfHits( aPatElectron.gsfTrack()->validFraction());

    float vertex0DistXY_;
    float vertex0DistZ_;

    if(recVtxs.isValid()) {

      int i_vertex = 0;
      for( reco::VertexCollection::const_iterator v_it=recVtxs->begin() ; v_it!=recVtxs->end() ; ++v_it ) {

	float distXY = aPatElectron.gsfTrack()->dxy(v_it->position());
	float distZ = aPatElectron.gsfTrack()->dz(v_it->position());
	float dist3D = sqrt(pow(distXY,2) + pow(distZ,2));

	if ( i_vertex == 0 ) {
	  vertex0DistXY_ = distXY;
	  vertex0DistZ_  = distZ ;
	  aElectron.setLeadVtxDistXY( vertex0DistXY_ );
	  aElectron.setLeadVtxDistZ( vertex0DistZ_ );
	}

	if( dist3D<minVtxDist3D ) {
	  minVtxDist3D = dist3D;
	  vertexIndex_ = int(std::distance(recVtxs->begin(),v_it));
	  vertexDistXY_ = distXY;
	  vertexDistZ_ = distZ;
	  aElectron.setVtxIndex( vertexIndex_ );
	  aElectron.setVtxDistXY( vertexDistXY_ );
	  aElectron.setVtxDistZ( vertexDistZ_ );

	}
	i_vertex++;
      }
    }

    aElectron.setscEta( aPatElectron.superCluster()->eta());
    aElectron.setscPhi( aPatElectron.superCluster()->phi());
    aElectron.setscPt( aPatElectron.superCluster()->energy() / cosh(aPatElectron.superCluster()->eta()));
    aElectron.setscRawEnergy( aPatElectron.superCluster()->rawEnergy());

    double dxy = aPatElectron.gsfTrack()->dxy(pv.position());
    aElectron.setdxy( dxy ) ;
    double dz = aPatElectron.gsfTrack()->dz(pv.position());
    aElectron.setdz( dz ) ;

    aElectron.setrho( rhoIso) ;
    
    aElectron.setPassConversionVeto( aPatElectron.passConversionVeto() );
    aElectron.setIsGsfCtfScPixChargeConsistent( aPatElectron.isGsfCtfScPixChargeConsistent());
    aElectron.setIsGsfScPixChargeConsistent( aPatElectron.isGsfScPixChargeConsistent());
    aElectron.setIsGsfCtfChargeConsistent( aPatElectron.isGsfCtfChargeConsistent());

    aElectron.setElectronIDs(aPatElectron.electronIDs());
    // for additional electron pids
    for (size_t i = 0; i < elecIDSrcs_.size(); ++i){
      ids[i].second = (*idhandles[i])[elecsRef];
      aElectron.setElectronID(ids[i]);
    }

    out->push_back(aElectron);
    ++j;
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

bool cat::CATElectronProducer::mcMatch( const reco::Candidate::LorentzVector& lepton, const edm::Handle<reco::GenParticleCollection> & genParticles )
{
  bool out = false;

  for (const reco::GenParticle & aGenPart : *genParticles){
    if( abs(aGenPart.pdgId()) != 11 ) continue;

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

bool cat::CATElectronProducer::MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact ) {
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
DEFINE_FWK_MODULE(CATElectronProducer);
