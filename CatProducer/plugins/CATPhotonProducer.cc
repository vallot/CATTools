#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "CATTools/DataFormats/interface/Photon.h"
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"



using namespace edm;
using namespace std;

namespace cat {

  class CATPhotonProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATPhotonProducer(const edm::ParameterSet & iConfig);

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
    int  mcMatch( const reco::Candidate::LorentzVector& lepton, const edm::Handle<reco::GenParticleCollection> & genParticles );
    double MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact );
    void findFirstNonPhotonMother(const reco::Candidate *particle,
				  int &ancestorPID, int &ancestorStatus);


    enum PhotonMatchType {UNMATCHED = 0,
                          MATCHED_FROM_GUDSCB,
                          MATCHED_FROM_PI0,
                          MATCHED_FROM_OTHER_SOURCES};

  private:

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<pat::Photon> > src_;
    edm::EDGetTokenT<edm::View<pat::Photon> > unsmearedPhotToken_;

    // Photon variables computed upstream in a special producer
    edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;
    
    edm::EDGetTokenT<double> rhoLabel_;

    bool runOnMC_;

    typedef std::pair<std::string, edm::InputTag> NameTag;
    typedef math::XYZPoint Point;
    std::vector<NameTag> phoIDSrcs_;
    std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > phoIDTokens_;
    const std::vector<std::string> photonIDs_;

  };

} // namespace

cat::CATPhotonProducer::CATPhotonProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("src"))),
  unsmearedPhotToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("unsmearedPhotons"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  mcLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"))),
  photonIDs_(iConfig.getParameter<std::vector<std::string> >("photonIDs"))
{
  produces<cat::PhotonCollection>();
  if (iConfig.existsAs<edm::ParameterSet>("photonIDSources")) {
    edm::ParameterSet idps = iConfig.getParameter<edm::ParameterSet>("photonIDSources");
    std::vector<std::string> names = idps.getParameterNamesForType<edm::InputTag>();
    for (std::vector<std::string>::const_iterator it = names.begin(), ed = names.end(); it != ed; ++it) {
      auto inputTag = idps.getParameter<edm::InputTag>(*it);
      phoIDSrcs_.push_back(NameTag(inputTag.instance(), inputTag));
    }
    phoIDTokens_ = edm::vector_transform(phoIDSrcs_, [this](NameTag const & tag){return mayConsume<edm::ValueMap<bool> >(tag.second);});
  }
}

void
cat::CATPhotonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  using namespace edm;
  using namespace std;

  runOnMC_ = !iEvent.isRealData();

  Handle<edm::View<pat::Photon> > src;
  iEvent.getByToken(src_, src);

  edm::Handle<edm::View<pat::Photon> > unsmearedPhotHandle;
  iEvent.getByToken(unsmearedPhotToken_, unsmearedPhotHandle);

  Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rhoIso = std::max(*(rhoHandle.product()), 0.0);

  Handle<reco::GenParticleCollection> genParticles;
  if (runOnMC_){
    iEvent.getByToken(mcLabel_,genParticles);
  }

  std::vector<edm::Handle<edm::ValueMap<bool> > > idhandles;
  std::vector<pat::Photon::IdPair>               ids;
  idhandles.resize(phoIDSrcs_.size());
  ids.resize(phoIDSrcs_.size());
  for (size_t i = 0; i < phoIDSrcs_.size(); ++i) {
    iEvent.getByToken(phoIDTokens_[i], idhandles[i]);
    ids[i].first = phoIDSrcs_[i].first;
  }
  

  std::unique_ptr<cat::PhotonCollection>  out(new cat::PhotonCollection());
  int j = 0;

  for (const pat::Photon & aPatPhoton : *src){
    cat::Photon aPhoton(aPatPhoton);
    auto phosRef = src->refAt(j);
    auto unsmearedPhotRef = unsmearedPhotHandle->refAt(j);

    if (runOnMC_){
      aPhoton.setGenParticleRef(aPatPhoton.genParticleRef());
      if(mcMatch(aPatPhoton.p4(), genParticles )  == 1) aPhoton.setMCMatched( true);
      else  aPhoton.setMCMatched( false);

      aPhoton.setSmearedScale(phosRef->pt()/unsmearedPhotRef->pt());
    }
    
    aPhoton.setRho(rhoIso);

    aPhoton.setPassElectronVeto(aPatPhoton.passElectronVeto() );
    aPhoton.setHasPixelSeed(aPatPhoton.hasPixelSeed() );
    // SC variables
    aPhoton.setSCrawEnergy(aPatPhoton.superCluster()->rawEnergy() );
    aPhoton.setSCPreShowerEnergy(aPatPhoton.superCluster()->preshowerEnergy() );
    aPhoton.setSCEta(aPatPhoton.superCluster()->eta() );
    aPhoton.setSCPhi(aPatPhoton.superCluster()->phi() );
    // add phi/eta width of SC?

    // shower variables + other
    aPhoton.setHoverE(aPatPhoton.hadTowOverEm()  );
    aPhoton.setR9(aPatPhoton.r9() );
    aPhoton.setSigmaiEtaiEta(aPatPhoton.full5x5_sigmaIetaIeta()) ;
    //aPhoton.setSigmaEtaEta(aPatPhoton.full5x5_sigmaEtaEta() );
    //aPhoton.setSigmaIphiIphi(aPatPhoton.full5x5_sigmaIphiIphi() ) ;

    //Energies
    //aPhoton.sete1_5(aPatPhoton.full5x5_e1x5() );
    //aPhoton.sete1_3(aPatPhoton.full5x5_e1x3() );
    //aPhoton.sete2_5Max(aPatPhoton.full5x5_e2x5Max() ); 
    //aPhoton.sete5_5(aPatPhoton.full5x5_e5x5() );
    //aPhoton.setfulr9(aPatPhoton.full5x5_r9() );
    
    aPhoton.setchargedHadronIso( aPatPhoton.chargedHadronIso() );
    aPhoton.setneutralHadronIso(aPatPhoton.neutralHadronIso() );
    aPhoton.setphotonIso(aPatPhoton.photonIso() );
    aPhoton.setpuhargedHadronIso(aPatPhoton.puChargedHadronIso() ); 
    
    if (phoIDSrcs_.size()){// for remade photon IDs
      for (size_t i = 0; i < phoIDSrcs_.size(); ++i){
        ids[i].second = (*idhandles[i])[unsmearedPhotRef];
        aPhoton.setPhotonID(ids[i]);
      }
    }
    else if (photonIDs_.size()){// for sphoted IDs in miniAOD
      for(unsigned int i = 0; i < photonIDs_.size(); i++){

        pat::Photon::IdPair pid(photonIDs_.at(i), aPatPhoton.photonID(photonIDs_.at(i)));
        aPhoton.setPhotonID(pid);
      }
    }
    else {
      aPhoton.setPhotonIDs(aPatPhoton.photonIDs());
    }
    
    aPhoton.setIsTight( aPhoton.photonID("tight") );
    aPhoton.setIsMedium( aPhoton.photonID("medium") );
    aPhoton.setIsLoose( aPhoton.photonID("loose") );
    aPhoton.setPassMVA( aPhoton.photonID("mva") );
      
    out->push_back(aPhoton);
    ++j;
  }

  iEvent.put(std::move(out));
}




int
cat::CATPhotonProducer::mcMatch( const reco::Candidate::LorentzVector& photon, const edm::Handle<reco::GenParticleCollection> & genParticles )
{

  double dR = 999;
  const reco::Candidate *closestPhoton = 0;


  for (const reco::GenParticle & aGenPart : *genParticles){
    const reco::Candidate *particle =  &aGenPart;
    // Drop everything that is not photon or not status 1
    if( abs(aGenPart.pdgId()) != 22 || aGenPart.status() != 1) continue;
    

    double dRtmp = MatchObjects(photon, aGenPart.p4(), false);
    
    if( dRtmp < dR ){
      dR = dRtmp;
      closestPhoton = particle;
    }
  }
  
  // See if the closest photon (if it exists) is close enough.
  // If not, no match found.
  if( !(closestPhoton != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }
  
  // Find ID of the parent of the found generator level photon match
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

  // Allowed parens: quarks pdgId 1-5, or a gluon 21
  std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
  if( !(std::find(allowedParents.begin(), 
		  allowedParents.end(), ancestorPID)
	!= allowedParents.end()) ){
    // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not. 
    if( abs(ancestorPID) == 111 )
      return MATCHED_FROM_PI0;
    else
      return MATCHED_FROM_OTHER_SOURCES;
  }
  return MATCHED_FROM_GUDSCB;
   
}

double
cat::CATPhotonProducer::MatchObjects( const reco::Candidate::LorentzVector& pasObj, const reco::Candidate::LorentzVector& proObj, bool exact )
{
  double proEta = proObj.eta();
  double proPhi = proObj.phi();
  double pasEta = pasObj.eta();
  double pasPhi = pasObj.phi();

  double dRval = deltaR(proEta, proPhi, pasEta, pasPhi);
  // If we are comparing two objects for which the candidates should
  // be exactly the same, cut hard. Otherwise take cuts from user.
  return  dRval;
}

void 
cat::CATPhotonProducer::findFirstNonPhotonMother(const reco::Candidate *particle,
						    int &ancestorPID, int &ancestorStatus){
  
  if( particle == 0 ){
    printf("PhotonNtuplerVIDDemo: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}


#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATPhotonProducer);
