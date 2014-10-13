/**
  \class    cat::CATJetProducer CATJetProducer.h "CATTools/CatProducer/interface/CATJetProducer.h"
  \brief    CAT Jet 
*/


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
      edm::EDGetTokenT<edm::View<pat::Jet> > src_;
      const std::vector<std::string> btagType_;

  };

} // namespace

cat::CATJetProducer::CATJetProducer(const edm::ParameterSet & iConfig) :
    src_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("src"))),
    btagType_(iConfig.getParameter<std::vector<std::string> >("btagType"))
{
    produces<std::vector<cat::Jet> >();
}

void 
cat::CATJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<pat::Jet> > src;
    iEvent.getByToken(src_, src);

    auto_ptr<vector<cat::Jet> >  out(new vector<cat::Jet>());

    for (View<pat::Jet>::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) {
      unsigned int idx = it - src->begin();
      const pat::Jet & aPatJet = src->at(idx);
 
      bool looseId = checkPFJetId( aPatJet ); 

      cat::Jet aJet(aPatJet);

      aJet.setLooseId( looseId );

      for(unsigned int i = 0; i < btagType_.size(); i++){
        const std::string tag(btagType_.at(i));
        aJet.setbTag( i, aPatJet.bDiscriminator(tag), tag);
      //  output->btag_[i] = input->bDiscriminator(tag);
      //  output->btagNames_[i] = tag;
      }
      
      //secondary vertex b-tagging information
      if( aPatJet.hasUserFloat("secvtxMass") ) aJet.setSecVtxMass( aPatJet.userFloat("secvtxMass") );
      if( aPatJet.hasUserFloat("Lxy") ) aJet.setLxy( aPatJet.userFloat("Lxy") );
      if( aPatJet.hasUserFloat("LxyErr") ) aJet.setLxyErr( aPatJet.userFloat("LxyErr") );

      Int_t partonFlavour = aPatJet.partonFlavour();
      aJet.setPartonFlavour( partonFlavour );
      Int_t partonPdgId = aPatJet.genParton() ? aPatJet.genParton()->pdgId() : 0;
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
