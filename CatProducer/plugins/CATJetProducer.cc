#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATJetProducer : public edm::EDProducer {
  public:
    explicit CATJetProducer(const edm::ParameterSet & iConfig);
    virtual ~CATJetProducer() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    bool checkPFJetId(const pat::Jet & jet);
      
    std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
    bool hasBottom(const reco::Candidate &c);
    bool hasCharm(const reco::Candidate &c);
    bool decayFromBHadron(const reco::Candidate &c);
    bool decayFromCHadron(const reco::Candidate &c);
    const reco::Candidate* lastBHadron(const reco::Candidate &c);
    const reco::Candidate* lastCHadron(const reco::Candidate &c);

  private:
    edm::EDGetTokenT<pat::JetCollection> src_;
    edm::EDGetTokenT<pat::JetCollection> shiftedEnDownSrc_;
    edm::EDGetTokenT<pat::JetCollection> shiftedEnUpSrc_;
    edm::EDGetTokenT<pat::JetCollection> smearedResSrc_;
    edm::EDGetTokenT<pat::JetCollection> smearedResDownSrc_;
    edm::EDGetTokenT<pat::JetCollection> smearedResUpSrc_;

    const std::vector<std::string> btagNames_;
    std::string uncertaintyTag_, payloadName_;
    bool runOnMC_;

  };

} // namespace

cat::CATJetProducer::CATJetProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  shiftedEnDownSrc_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("shiftedEnDownSrc"))),
  shiftedEnUpSrc_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("shiftedEnUpSrc"))),
  smearedResSrc_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("smearedResSrc"))),
  smearedResDownSrc_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("smearedResDownSrc"))),
  smearedResUpSrc_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("smearedResUpSrc"))),
  btagNames_(iConfig.getParameter<std::vector<std::string> >("btagNames"))
{
  produces<std::vector<cat::Jet> >();
  //uncertaintyTag_    = iConfig.getParameter<std::string>("UncertaintyTag");
}

void 
cat::CATJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  runOnMC_ = !iEvent.isRealData();

  edm::Handle<pat::JetCollection> src;
  iEvent.getByToken(src_, src);

  edm::Handle<pat::JetCollection> shiftedEnDownSrc;
  edm::Handle<pat::JetCollection> shiftedEnUpSrc;
  edm::Handle<pat::JetCollection> smearedResSrc;
  edm::Handle<pat::JetCollection> smearedResDownSrc;
  edm::Handle<pat::JetCollection> smearedResUpSrc;
  if (runOnMC_){
    iEvent.getByToken(shiftedEnDownSrc_, shiftedEnDownSrc);
    iEvent.getByToken(shiftedEnUpSrc_, shiftedEnUpSrc);
    iEvent.getByToken(smearedResSrc_, smearedResSrc);
    iEvent.getByToken(smearedResDownSrc_, smearedResDownSrc);
    iEvent.getByToken(smearedResUpSrc_, smearedResUpSrc);
  }
  
  auto_ptr<vector<cat::Jet> >  out(new vector<cat::Jet>());
  int j = 0;
  for (const pat::Jet &aPatJet : *src) {
    bool looseId = checkPFJetId( aPatJet );
    cat::Jet aJet(aPatJet);

    if (runOnMC_){
      aJet.setShiftedEnDown(shiftedEnDownSrc->at(j).pt() );
      aJet.setShiftedEnUp(shiftedEnUpSrc->at(j).pt() );
      aJet.setSmearedRes(smearedResSrc->at(j).pt() );
      aJet.setSmearedResDown(smearedResDownSrc->at(j).pt() );
      aJet.setSmearedResUp(smearedResUpSrc->at(j).pt() );
      // adding genJet
      aJet.setGenJetRef(aPatJet.genJetFwdRef());
    }
    ++j;
    aJet.setLooseId( looseId );
    if( aPatJet.hasUserFloat("pileupJetId:fullDiscriminant") )
      aJet.setPileupJetId( aPatJet.userFloat("pileupJetId:fullDiscriminant") );

    if (btagNames_.size() == 0){
      aJet.setBDiscriminators(aPatJet.getPairDiscri());
    }
    else {
      for(unsigned int i = 0; i < btagNames_.size(); i++){
	aJet.addBDiscriminatorPair(std::make_pair(btagNames_.at(i), aPatJet.bDiscriminator(btagNames_.at(i)) ));
      }
    }

    //secondary vertex b-tagging information
    if( aPatJet.hasUserFloat("vtxMass") ) aJet.setVtxMass( aPatJet.userFloat("vtxMass") );
    if( aPatJet.hasUserFloat("vtxNtracks") ) aJet.setVtxNtracks( aPatJet.userFloat("vtxNtracks") );
    if( aPatJet.hasUserFloat("vtx3DVal") ) aJet.setVtx3DVal( aPatJet.userFloat("vtx3DVal") );
    if( aPatJet.hasUserFloat("vtx3DSig") ) aJet.setVtx3DSig( aPatJet.userFloat("vtx3DSig") );

    aJet.setHadronFlavour(aPatJet.hadronFlavour());
    aJet.setPartonFlavour(aPatJet.partonFlavour());
    int partonPdgId = aPatJet.genParton() ? aPatJet.genParton()->pdgId() : 0;
    aJet.setPartonPdgId(partonPdgId);

    out->push_back(aJet);
  }

  iEvent.put(out);
}

bool cat::CATJetProducer::checkPFJetId(const pat::Jet & jet){
  //Loose PF Jet id
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  //debug
  bool out = false;
  if( (jet.neutralHadronEnergy() + jet.HFHadronEnergy() ) / jet.energy() < 0.99
      &&jet.neutralEmEnergyFraction() < 0.99
      &&jet.numberOfDaughters() > 1
      &&(jet.chargedHadronEnergyFraction() > 0 || abs(jet.eta()) > 2.4)
      &&(jet.chargedMultiplicity() > 0 || abs(jet.eta()) > 2.4)
      &&(jet.chargedEmEnergyFraction() < 0.99 || abs(jet.eta()) > 2.4)
      ) out = true;

  return out;
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATJetProducer);
