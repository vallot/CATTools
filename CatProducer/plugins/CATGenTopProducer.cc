#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace edm;
using namespace std;

namespace cat {

  class CATGenTopProducer : public edm::stream::EDProducer<> {
  public:
    explicit CATGenTopProducer(const edm::ParameterSet & iConfig);

    void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

  private:
    edm::EDGetTokenT<reco::GenJetCollection> genJetLabel_;
    edm::EDGetTokenT<reco::GenParticleCollection> mcParticleLabel_;

    const edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken_;
    const edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken_;

    const edm::EDGetTokenT<std::vector<int> > genCHadJetIndexToken_;
    const edm::EDGetTokenT<std::vector<int> > genCHadFlavourToken_;


  };

} // namespace

cat::CATGenTopProducer::CATGenTopProducer(const edm::ParameterSet & iConfig) :
  genJetLabel_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetLabel"))),
  mcParticleLabel_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcParticleLabel"))),
  genBHadJetIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadJetIndex"))),
  genBHadFlavourToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadFlavour"))),
  genCHadJetIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadJetIndex"))),
  genCHadFlavourToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadFlavour")))
{
  produces<std::vector<cat::GenTop> >();
}

void
cat::CATGenTopProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetLabel_, genJets);

  Handle<reco::GenParticleCollection> mcParticles;
  iEvent.getByToken(mcParticleLabel_, mcParticles);

  // Access B hadrons information
  Handle<std::vector<int> > genBHadFlavour;
  iEvent.getByToken(genBHadFlavourToken_, genBHadFlavour);

  Handle<std::vector<int> > genBHadJetIndex;
  iEvent.getByToken(genBHadJetIndexToken_, genBHadJetIndex);

  // Access C hadrons information
  Handle<std::vector<int> > genCHadFlavour;
  iEvent.getByToken(genCHadFlavourToken_, genCHadFlavour);

  Handle<std::vector<int> > genCHadJetIndex;
  iEvent.getByToken(genCHadJetIndexToken_, genCHadJetIndex);

  cat::GenTop aGenTop;
  aGenTop.building( genJets, mcParticles, genBHadFlavour, genBHadJetIndex, genCHadFlavour, genCHadJetIndex );

  unique_ptr<vector<cat::GenTop> >  out(new vector<cat::GenTop>());

  out->push_back(aGenTop);

  iEvent.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATGenTopProducer);
