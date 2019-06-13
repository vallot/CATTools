#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
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
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "FWCore/Utilities/interface/transform.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"


#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "CATTools/CommonTools/interface/GenParticleHelper.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATElectronProducer : public edm::stream::EDProducer<>
  {
  public:
    explicit CATElectronProducer(const edm::ParameterSet & iConfig);

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
    bool mcMatch( const reco::Candidate::LorentzVector& lepton, const std::vector<reco::GenParticleRef>& genParticles ) const;
    double getMiniRelIso(edm::Handle<pat::PackedCandidateCollection> pfcands,  const reco::Candidate::LorentzVector& ptcl, double  r_iso_min, double r_iso_max, double kt_scale,double rhoIso, double AEff) const;

  private:

    float getEffArea( float dR, float scEta );
    edm::EDGetTokenT<edm::View<pat::Electron> > src_;
    edm::EDGetTokenT<edm::View<pat::Electron> > unsmearedElecToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
    edm::EDGetTokenT<pat::PackedCandidateCollection>        pfSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamLineSrc_;
    edm::EDGetTokenT<double> rhoLabel_;

    bool runOnMC_;

    typedef std::pair<std::string, edm::InputTag> NameTag;
    typedef math::XYZPoint Point;

    std::vector<NameTag> elecIDSrcs_;
    std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > elecIDTokens_;
    const std::vector<std::string> electronIDs_;

  };

} // namespace



double cat::CATElectronProducer::getMiniRelIso(edm::Handle<pat::PackedCandidateCollection> pfcands,
					       const reco::Candidate::LorentzVector& ptcl,
					       double r_iso_min, double r_iso_max, double kt_scale, double rhoIso, double AEff) const
{

  if (ptcl.pt()<5.) return 99999.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if (fabs(ptcl.eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  const double ptThresh(0.);
  const double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl.pt()));

  double iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
  for (const pat::PackedCandidate &pfc : *pfcands) {
    const unsigned int absId = std::abs(pfc.pdgId());
    if ( absId<7 ) continue;

    const double dr = deltaR(pfc, ptcl);
    if (dr > r_iso) continue;

    if (pfc.charge()==0){
      if (pfc.pt() <= ptThresh) continue;

      if (abs(pfc.pdgId())==22) {
        if(dr >= deadcone_ph) iso_ph += pfc.pt();
      }
      else if (abs(pfc.pdgId())==130) {
        if(dr >= deadcone_nh) iso_nh += pfc.pt();
      }
    }
    else if (pfc.fromPV()>1){
      if (abs(pfc.pdgId())==211 and dr >= deadcone_ch ) iso_ch += pfc.pt();
    }
    else {
      if (pfc.pt()>ptThresh and dr >= deadcone_pu) iso_pu += pfc.pt();
    }
  }

  const double conesize_correction = r_iso*r_iso/0.09;
  const double iso = iso_ch  + std::max(0.0, iso_nh + iso_ph - rhoIso*AEff*conesize_correction);
  return iso/ptcl.pt();
}

cat::CATElectronProducer::CATElectronProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("src"))),
  unsmearedElecToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("unsmaredElectrons"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  mcLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  pfSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  beamLineSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamLineSrc"))),
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"))),
  electronIDs_(iConfig.getParameter<std::vector<std::string> >("electronIDs"))
{
  produces<cat::ElectronCollection>();

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

void cat::CATElectronProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  using namespace edm;
  using namespace std;

  runOnMC_ = !iEvent.isRealData();

  edm::Handle<edm::View<pat::Electron> > src;
  iEvent.getByToken(src_, src);

  edm::Handle<edm::View<pat::Electron> > unsmearedElecHandle;
  iEvent.getByToken(unsmearedElecToken_, unsmearedElecHandle);

  edm::Handle<reco::GenParticleCollection> genParticles;
  std::vector<reco::GenParticleRef> genElectrons;
  if (runOnMC_){
    iEvent.getByToken(mcLabel_,genParticles);
    genElectrons = cat::selectGenParticles(genParticles, 11, {23,24});
  }

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_, recVtxs);
  reco::Vertex pv;
  if (recVtxs->size()) pv = recVtxs->at(0);
  const GlobalPoint pVertex(pv.position().x(),pv.position().y(),pv.position().z());

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfSrc_, pfcands);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rhoIso = std::max(*(rhoHandle.product()), 0.0);

  std::vector<edm::Handle<edm::ValueMap<bool> > > idhandles;
  std::vector<pat::Electron::IdPair>               ids;
  idhandles.resize(elecIDSrcs_.size());
  ids.resize(elecIDSrcs_.size());
  for (size_t i = 0; i < elecIDSrcs_.size(); ++i) {
    iEvent.getByToken(elecIDTokens_[i], idhandles[i]);
    ids[i].first = elecIDSrcs_[i].first;
  }

  std::unique_ptr<cat::ElectronCollection>  out(new cat::ElectronCollection());
  for ( int j=0, n=src->size(); j<n; ++j ) {
    const auto aPatElectron = src->at(j);
    cat::Electron aElectron(aPatElectron);
    auto elecsRef = src->refAt(j);
    auto unsmearedElecRef = unsmearedElecHandle->refAt(j);
    // nan protection - smearing fails for soft electrons
    if ( std::isnan(std::abs(aElectron.p())) ) aElectron = *unsmearedElecRef;

    if (runOnMC_){
      aElectron.setGenParticleRef(aPatElectron.genParticleRef());
      aElectron.setMCMatched( mcMatch( aElectron.p4(), genElectrons ) );
    }
    aElectron.setSmearedScale( aPatElectron.userFloat("ecalTrkEnergyPostCorr") / aPatElectron.energy());
    aElectron.setSmearedScaleUnc( aPatElectron.userFloat("energyScaleUp"),
                                  aPatElectron.userFloat("energyScaleDown"),
                                  aPatElectron.userFloat("energySigmaUp"),
                                  aPatElectron.userFloat("energySigmaDown"));
    aElectron.setIsGsfCtfScPixChargeConsistent( aPatElectron.isGsfCtfScPixChargeConsistent() );
    aElectron.setIsEB( aPatElectron.isEB() );

    aElectron.setChargedHadronIso04( aPatElectron.chargedHadronIso() );
    aElectron.setNeutralHadronIso04( aPatElectron.neutralHadronIso() );
    aElectron.setPhotonIso04( aPatElectron.photonIso() );
    aElectron.setPUChargedHadronIso04( aPatElectron.puChargedHadronIso() );

    aElectron.setChargedHadronIso03( aPatElectron.pfIsolationVariables().sumChargedHadronPt );
    aElectron.setNeutralHadronIso03( aPatElectron.pfIsolationVariables().sumNeutralHadronEt );
    aElectron.setPhotonIso03( aPatElectron.pfIsolationVariables().sumPhotonEt );
    aElectron.setPUChargedHadronIso03( aPatElectron.pfIsolationVariables().sumPUPt );

    const float scEta = aPatElectron.superCluster()->eta();
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
    aElectron.setMiniRelIso(getMiniRelIso( pfcands, aElectron.p4(), 0.05, 0.2, 10., rhoIso,elEffArea03));


    aElectron.setscEta(scEta);
    aElectron.setPassConversionVeto( aPatElectron.passConversionVeto() );

    if (elecIDSrcs_.size()){// for remade electron IDs
      for (size_t i = 0; i < elecIDSrcs_.size(); ++i){
        ids[i].second = (*idhandles[i])[unsmearedElecRef];
        aElectron.setElectronID(ids[i]);
      }
    }
    else if (electronIDs_.size()){// for selected IDs in miniAOD
      for(unsigned int i = 0; i < electronIDs_.size(); i++){
        pat::Electron::IdPair pid(electronIDs_.at(i), aPatElectron.electronID(electronIDs_.at(i)));
        aElectron.setElectronID(pid);
      }
    }
    else {
      aElectron.setElectronIDs(aPatElectron.electronIDs());
    }

    aElectron.setIsPF( aPatElectron.isPF() );
    aElectron.setIsTight( aElectron.electronID("cutBasedElectronID-Fall17-94X-V2-tight") );
    aElectron.setIsMedium( aElectron.electronID("cutBasedElectronID-Fall17-94X-V2-medium") );
    aElectron.setIsLoose( aElectron.electronID("cutBasedElectronID-Fall17-94X-V2-loose") );

    reco::GsfTrackRef theTrack = aPatElectron.gsfTrack();
    aElectron.setDxy( theTrack->dxy(pv.position()) );
    aElectron.setDz( theTrack->dz(pv.position()) );
    aElectron.setVertex(Point(theTrack->vx(),theTrack->vy(),theTrack->vz()));

    reco::TransientTrack eletranstrack = trackBuilder->build(theTrack);
    TrajectoryStateOnSurface eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
    if (eleTSOS.isValid()) {
      std::pair<bool, Measurement1D>     eleIPpair = IPTools::signedTransverseImpactParameter(eletranstrack, eleTSOS.globalDirection(), pv);

      float eleSignificanceIP = eleIPpair.second.significance();
      aElectron.setIpSignficance(eleSignificanceIP);
    }
    out->push_back(aElectron);
  }
  iEvent.put(std::move(out));
}

float
cat::CATElectronProducer::getEffArea( float dR, float scEta)
{
  // https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
  float absEta = std::abs(scEta);
  if ( 0.0000 >= absEta && absEta < 1.0000 ) return 0.1440;
  if ( 1.0000 >= absEta && absEta < 1.4790 ) return 0.1562;
  if ( 1.4790 >= absEta && absEta < 2.0000 ) return 0.1032;
  if ( 2.0000 >= absEta && absEta < 2.2000 ) return 0.0859;
  if ( 2.2000 >= absEta && absEta < 2.3000 ) return 0.1116;
  if ( 2.3000 >= absEta && absEta < 2.4000 ) return 0.1321;
  if ( 2.4000 >= absEta && absEta < 5.0000 ) return 0.1654;
  return 0;
}

bool cat::CATElectronProducer::mcMatch( const reco::Candidate::LorentzVector& lepton, const std::vector<reco::GenParticleRef>& genElectrons) const
{
  for (const auto& aGenPart : genElectrons){
    if ( cat::isMatchedByDRDPt(lepton, aGenPart->p4(), false) ) return true;
  }
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATElectronProducer);
