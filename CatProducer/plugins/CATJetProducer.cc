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
    bool checkPFJetIdTight(const pat::Jet & jet);
    void getJER(const double jetEta, double& cJER, double& cJERUp, double& cJERDn) const;
      
    std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
    bool hasBottom(const reco::Candidate &c);
    bool hasCharm(const reco::Candidate &c);
    bool decayFromBHadron(const reco::Candidate &c);
    bool decayFromCHadron(const reco::Candidate &c);
    const reco::Candidate* lastBHadron(const reco::Candidate &c);
    const reco::Candidate* lastCHadron(const reco::Candidate &c);

  private:
    edm::EDGetTokenT<pat::JetCollection> src_;

    const std::vector<std::string> btagNames_;
    std::string uncertaintyTag_, payloadName_;
    bool runOnMC_;
    //PFJetIDSelectionFunctor pfjetIDFunctor;
    JetCorrectionUncertainty *jecUnc;
  };

} // namespace

cat::CATJetProducer::CATJetProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  btagNames_(iConfig.getParameter<std::vector<std::string> >("btagNames")),
  payloadName_(iConfig.getParameter<std::string>("payloadName"))
{
  produces<std::vector<cat::Jet> >();
  ///  pfjetIDFunctor = PFJetIDSelectionFunctor(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::LOOSE);
}

void 
cat::CATJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  runOnMC_ = !iEvent.isRealData();

  edm::Handle<pat::JetCollection> src;
  iEvent.getByToken(src_, src);

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  if (payloadName_.size()){
    // temp measure - payloadName should be AK4PFchs, but PHYS14_25_V2 does not have uncertainty 
    iSetup.get<JetCorrectionsRecord>().get(payloadName_,JetCorParColl); 
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }

  auto_ptr<vector<cat::Jet> >  out(new vector<cat::Jet>());
  for (const pat::Jet &aPatJet : *src) {

    cat::Jet aJet(aPatJet);

    aJet.addJecFactorPair(std::make_pair("Uncorrected", aPatJet.jecFactor("Uncorrected") ) );
    aJet.addJecFactorPair(std::make_pair("L1FastJet", aPatJet.jecFactor("L1FastJet") ) );
    aJet.addJecFactorPair(std::make_pair("L2Relative", aPatJet.jecFactor("L2Relative") ) );
    aJet.addJecFactorPair(std::make_pair("L3Absolute", aPatJet.jecFactor("L3Absolute") ) );
    // aJet.addJecFactorPair(std::make_pair("L5Flavor_gJ", aPatJet.jecFactor("L5Flavor_gJ") ) );
    // aJet.addJecFactorPair(std::make_pair("L5Flavor_qT", aPatJet.jecFactor("L5Flavor_qT") ) );
    // aJet.addJecFactorPair(std::make_pair("L5Flavor_cT", aPatJet.jecFactor("L5Flavor_cT") ) );
    // aJet.addJecFactorPair(std::make_pair("L5Flavor_bT", aPatJet.jecFactor("L5Flavor_bT") ) );

    if (aPatJet.isPFJet() )
      aJet.setPFSpecific( aPatJet.pfSpecific() );
    
    bool looseId = checkPFJetId( aPatJet );
    bool tightId = checkPFJetIdTight( aPatJet );
    aJet.setLooseId( looseId );
    aJet.setTightId( tightId );

    if( aPatJet.hasUserFloat("pileupJetId:fullDiscriminant") )
      aJet.setPileupJetId( aPatJet.userFloat("pileupJetId:fullDiscriminant") );

    if (btagNames_.size() == 0){
      aJet.setBDiscriminators(aPatJet.getPairDiscri());
      // const std::vector<std::pair<std::string, float> > bpair = aPatJet.getPairDiscri();
      // for (unsigned int bb =0; bb < bpair.size(); bb++){cout << bpair[bb].first <<endl;}
    }
    else {
      for(unsigned int i = 0; i < btagNames_.size(); i++){
	aJet.addBDiscriminatorPair(std::make_pair(btagNames_.at(i), aPatJet.bDiscriminator(btagNames_.at(i)) ));
      }
    }
    //cout << "jet pt " << aJet.pt() <<" eta " << aJet.eta() <<endl;

    //secondary vertex b-tagging information
    if( aPatJet.hasUserFloat("vtxMass") ) aJet.setVtxMass( aPatJet.userFloat("vtxMass") );
    if( aPatJet.hasUserFloat("vtxNtracks") ) aJet.setVtxNtracks( aPatJet.userFloat("vtxNtracks") );
    if( aPatJet.hasUserFloat("vtx3DVal") ) aJet.setVtx3DVal( aPatJet.userFloat("vtx3DVal") );
    if( aPatJet.hasUserFloat("vtx3DSig") ) aJet.setVtx3DSig( aPatJet.userFloat("vtx3DSig") );

    aJet.setHadronFlavour(aPatJet.hadronFlavour());
    aJet.setPartonFlavour(aPatJet.partonFlavour());
    int partonPdgId = aPatJet.genParton() ? aPatJet.genParton()->pdgId() : 0;
    aJet.setPartonPdgId(partonPdgId);

    // setting JEC uncertainty
    if (payloadName_.size()){
      jecUnc->setJetEta(aJet.eta());
      jecUnc->setJetPt(aJet.pt()); // here you must use the CORRECTED jet pt
      double unc = jecUnc->getUncertainty(true);
      aJet.setShiftedEnUp( (1. + unc) );
      jecUnc->setJetEta(aJet.eta());
      jecUnc->setJetPt(aJet.pt()); // here you must use the CORRECTED jet pt
      unc = jecUnc->getUncertainty(false);
      aJet.setShiftedEnDown( (1. - unc) );
    }
    float fJER   = 0.;
    float fJERUp = 0.;
    float fJERDn = 0.;
    if (runOnMC_){
      // adding genJet
      aJet.setGenJetRef(aPatJet.genJetFwdRef());
      aJet.setGenParticleRef(aPatJet.genParticleRef());

      // setting JES 
      if ( aPatJet.genJet() ){
	double cJER, cJERUp, cJERDn;
	getJER(aJet.eta(), cJER, cJERUp, cJERDn);

	const double jetPt = aJet.pt();
	const double genJetPt = aPatJet.genJet()->pt();
	const double dPt = jetPt-genJetPt;

	fJER   = max(0., (genJetPt+dPt*cJER  )/jetPt);
	fJERUp = max(0., (genJetPt+dPt*cJERUp)/jetPt);
	fJERDn = max(0., (genJetPt+dPt*cJERDn)/jetPt);
      }
    }
    aJet.setSmearedRes(fJER);
    aJet.setSmearedResDown(fJERDn);
    aJet.setSmearedResUp(fJERUp);

    out->push_back(aJet);
  }

  iEvent.put(out);
}

bool cat::CATJetProducer::checkPFJetId(const pat::Jet & jet){
  //Loose PF Jet id
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  //debug
  bool out = false;
  if ( jet.neutralHadronEnergyFraction() < 0.99
       &&jet.neutralEmEnergyFraction() < 0.99
       &&jet.numberOfDaughters() > 1
       &&(jet.chargedHadronEnergyFraction() > 0 || abs(jet.eta()) > 2.4)
       &&(jet.chargedMultiplicity() > 0 || abs(jet.eta()) > 2.4)
       &&(jet.chargedEmEnergyFraction() < 0.99 || abs(jet.eta()) > 2.4)
       ) out = true;

  return out;
}

bool cat::CATJetProducer::checkPFJetIdTight(const pat::Jet & jet){
  //Tight PF Jet id
  bool out = false;
  if ( jet.neutralHadronEnergyFraction() < 0.90
       &&jet.neutralEmEnergyFraction() < 0.90
       &&jet.numberOfDaughters() > 1
       &&(jet.chargedHadronEnergyFraction() > 0 || abs(jet.eta()) > 2.4)
       &&(jet.chargedMultiplicity() > 0 || abs(jet.eta()) > 2.4)
       &&(jet.chargedEmEnergyFraction() < 0.90 || abs(jet.eta()) > 2.4)
       ) out = true;
  return out;
}

void cat::CATJetProducer::getJER(const double jetEta, double& cJER, double& cJERUp, double& cJERDn) const{
  // 2012 values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
  // need to update for run2, must be better way to impliment these corrections
  const double absEta = std::abs(jetEta);
  if      ( absEta < 0.5 ) { cJER = 1.079; cJERUp = 1.105; cJERDn = 1.053; }
  else if ( absEta < 1.1 ) { cJER = 1.099; cJERUp = 1.127; cJERDn = 1.071; }
  else if ( absEta < 1.7 ) { cJER = 1.121; cJERUp = 1.150; cJERDn = 1.092; }
  else if ( absEta < 2.3 ) { cJER = 1.208; cJERUp = 1.254; cJERDn = 1.162; }
  else if ( absEta < 2.8 ) { cJER = 1.254; cJERUp = 1.316; cJERDn = 1.192; }
  else if ( absEta < 3.2 ) { cJER = 1.395; cJERUp = 1.458; cJERDn = 1.332; }
  else if ( absEta < 5.0 ) { cJER = 1.056; cJERUp = 1.247; cJERDn = 0.865; }
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATJetProducer);
