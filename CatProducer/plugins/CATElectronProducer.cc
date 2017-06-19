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
    int getSNUID(float, float, float, float, float, float, int, bool, float);
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
      aElectron.setSmearedScale(aElectron.pt()/unsmearedElecRef->pt());
    }
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
    aElectron.setIsTight( aElectron.electronID("tight") );
    aElectron.setIsMedium( aElectron.electronID("medium") );
    aElectron.setIsLoose( aElectron.electronID("loose") );

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


    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
    // The if protects against ecalEnergy == inf or zero
    // (always the case for miniAOD for electrons <5 GeV)
    float eoverp = 1e30;
    if( aPatElectron.ecalEnergy() != 0 and std::isfinite(aPatElectron.ecalEnergy()) ) {
      eoverp = std::abs(1.0/aPatElectron.ecalEnergy() - aPatElectron.eSuperClusterOverP()/aPatElectron.ecalEnergy() ) ;
    }

    int snu_id = getSNUID(aPatElectron.full5x5_sigmaIetaIeta(), abs(aPatElectron.deltaEtaSuperClusterTrackAtVtx() ), abs(aPatElectron.deltaPhiSuperClusterTrackAtVtx() ), aPatElectron.hcalOverEcal(), eoverp, abs(aElectron.dz()) , aPatElectron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS), aPatElectron.passConversionVeto(),scEta);
    aElectron.setSNUID(snu_id);

    // Fill the validity flag of triggered MVA
    bool isTrigMVAValid = false;
    const double pt = aElectron.pt();
    if ( pt > 15 ) {
      const double abseta = std::abs(scEta);
      const double full5x5_sigmaIetaIeta = aPatElectron.full5x5_sigmaIetaIeta();
      const double hcalOverEcal =  aPatElectron.hcalOverEcal();
      const double ecalPFClusterIso = aPatElectron.ecalPFClusterIso();
      const double hcalPFClusterIso = aPatElectron.hcalPFClusterIso();
      const double dr03TkSumPt = aPatElectron.dr03TkSumPt();
      if ( abseta < 1.479 ) { // Barrel
        const double deltaEtaSuperClusterTrackAtVtx = aPatElectron.deltaEtaSuperClusterTrackAtVtx();
        const double deltaPhiSuperClusterTrackAtVtx = aPatElectron.deltaPhiSuperClusterTrackAtVtx();
        if ( full5x5_sigmaIetaIeta < 0.012 and
             hcalOverEcal < 0.09 and
             ecalPFClusterIso < 0.37*pt and
             hcalPFClusterIso < 0.25*pt and
             dr03TkSumPt < 0.18*pt and
             std::abs(deltaEtaSuperClusterTrackAtVtx) < 0.0095 and
             std::abs(deltaPhiSuperClusterTrackAtVtx) < 0.065 ) isTrigMVAValid = true;
      }
      else { // Endcap
        if ( full5x5_sigmaIetaIeta < 0.033 and
             hcalOverEcal < 0.09 and
             ecalPFClusterIso < 0.45*pt and
             hcalPFClusterIso < 0.28*pt and
             dr03TkSumPt < 0.18*pt ) isTrigMVAValid = true;
      }
    }
    aElectron.setTrigMVAValid(isTrigMVAValid);

    out->push_back(aElectron);
  }
  iEvent.put(std::move(out));
}


int cat::CATElectronProducer::getSNUID(float full5x5_sigmaIetaIeta, float deltaEtaSuperClusterTrackAtVtx, float deltaPhiSuperClusterTrackAtVtx, float hoverE, float eoverp, float dz, int exp_miss_innerhits, bool pass_conversion_veto, float sceta){

  //----------------------------------------------------------------------
  // Barrel electron cut values
  //----------------------------------------------------------------------
  //Spring15 selection, 25ns selection
  //string id [4] = {"veto", "loose","medium", "tight" };

  double l_b_sieie   [4] = { 0.0115, 0.011, 0.00998, 0.00998};
  double l_b_dEtaIn  [4] = { 0.00749,0.00477, 0.00311, 0.00308};
  double l_b_dPhiIn  [4] = { 0.228, 0.222, 0.103, 0.0816};
  double l_b_hoe     [4] = { 0.356, 0.298, 0.253, 0.0414};
  double l_b_ep      [4] = { 0.299, 0.241, 0.134, 0.0129};
  int    l_b_missHits[4] = { 2, 1, 1, 1};

  //----------------------------------------------------------------------
  // Endcap electron cut values
  //----------------------------------------------------------------------

  double l_e_sieie   [4] = { 0.037, 0.0314, 0.0298, 0.0292};
  double l_e_dEtaIn  [4] = { 0.00895, 0.00868, 0.00609, 0.00605};
  double l_e_dPhiIn  [4] = { 0.213, 0.213, 0.045, 0.0394};
  double l_e_hoe     [4] = { 0.211, 0.101, 0.0878, 0.0641};
  double l_e_ep      [4] = { 0.15, 0.14, 0.13,0.0129};
  int    l_e_missHits[4] = { 3, 1, 1, 1};

  int flag_id=0;
  for(int i=0; i < 4; i++){
    bool pass_id=true;
    if ( std::abs(sceta) < 1.479 ){
      if(full5x5_sigmaIetaIeta >= l_b_sieie[i])pass_id = false;
      if(deltaEtaSuperClusterTrackAtVtx >= l_b_dEtaIn[i])pass_id = false;
      if(deltaPhiSuperClusterTrackAtVtx >= l_b_dPhiIn[i])pass_id = false;
      if(hoverE >= l_b_hoe[i])pass_id = false;
      if(eoverp >= l_b_ep[i])pass_id = false;
      if(exp_miss_innerhits > l_b_missHits[i])pass_id = false;
      if(!pass_conversion_veto) pass_id = false;
    }
    else   if ( std::abs(sceta) < 2.5 ){
      if(full5x5_sigmaIetaIeta >= l_e_sieie[i])pass_id = false;
      if(deltaEtaSuperClusterTrackAtVtx >= l_e_dEtaIn[i])pass_id = false;
      if(deltaPhiSuperClusterTrackAtVtx>= l_e_dPhiIn[i])pass_id = false;
      if(hoverE>= l_e_hoe[i])pass_id = false;
      if(eoverp>= l_e_ep[i])pass_id = false;
      if(exp_miss_innerhits > l_e_missHits[i])pass_id = false;
      if(!pass_conversion_veto) pass_id = false;
    }
    if(pass_id)flag_id += pow(10,i);
  }

  return flag_id;
}

float
cat::CATElectronProducer::getEffArea( float dR, float scEta)
{
  // ElectronEffectiveArea::ElectronEffectiveAreaTarget electronEATarget;
  // if ( runOnMC_ ) electronEATarget = ElectronEffectiveArea::kEleEAFall11MC;
  // else electronEATarget = ElectronEffectiveArea::kEleEAData2012;
  // if( dR < 0.35)
  //   return ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, scEta, electronEATarget);
  // else
  //   return ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, scEta, electronEATarget);

  // new effArea  https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
  float absEta = std::abs(scEta);
  if ( 0.0000 >= absEta && absEta < 1.0000 ) return 0.1703;
  if ( 1.0000 >= absEta && absEta < 1.4790 ) return 0.1715;
  if ( 1.4790 >= absEta && absEta < 2.0000 ) return 0.1213;
  if ( 2.0000 >= absEta && absEta < 2.2000 ) return 0.1230;
  if ( 2.2000 >= absEta && absEta < 2.3000 ) return 0.1635;
  if ( 2.3000 >= absEta && absEta < 2.4000 ) return 0.1937;
  if ( 2.4000 >= absEta && absEta < 5.0000 ) return 0.2393;
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
