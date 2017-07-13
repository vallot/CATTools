#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"

using namespace edm;
using namespace std;

namespace cat {

class CATJetProducer : public edm::stream::EDProducer<>
{
public:
  explicit CATJetProducer(const edm::ParameterSet & iConfig);

  void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;

  std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
  bool hasBottom(const reco::Candidate &c);
  bool hasCharm(const reco::Candidate &c);
  bool decayFromBHadron(const reco::Candidate &c);
  bool decayFromCHadron(const reco::Candidate &c);
  const reco::Candidate* lastBHadron(const reco::Candidate &c);
  const reco::Candidate* lastCHadron(const reco::Candidate &c);

private:
  edm::EDGetTokenT<pat::JetCollection> src_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> qgToken_;
  std::vector<std::string> flavTagNames_;
  std::vector<edm::EDGetTokenT<edm::ValueMap<float>>> flavTagTokens_;

  const std::vector<std::string> btagNames_;
  std::string uncertaintyTag_, payloadName_;
  const std::string jetResFilePath_, jetResSFFilePath_;
  bool setGenParticle_;
  bool runOnMC_;
  //PFJetIDSelectionFunctor pfjetIDFunctor;
  JetCorrectionUncertainty *jecUnc;

  CLHEP::HepRandomEngine* rng_;
};

} // namespace

cat::CATJetProducer::CATJetProducer(const edm::ParameterSet & iConfig) :
  src_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  qgToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("qgLikelihood"))),
  btagNames_(iConfig.getParameter<std::vector<std::string> >("btagNames")),
  payloadName_(iConfig.getParameter<std::string>("payloadName")),
  jetResFilePath_(edm::FileInPath(iConfig.getParameter<std::string>("jetResFile")).fullPath()),
  jetResSFFilePath_(edm::FileInPath(iConfig.getParameter<std::string>("jetResSFFile")).fullPath()),
  setGenParticle_(iConfig.getParameter<bool>("setGenParticle"))
{
  for ( auto label : iConfig.getParameter<std::vector<edm::InputTag>>("flavTagLabels") ) {
    const std::string name = label.label() + ":" + label.instance();
    flavTagNames_.push_back(name);
    flavTagTokens_.push_back(consumes<edm::ValueMap<float>>(label));
  }

  produces<cat::JetCollection>();
}

void cat::CATJetProducer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&)
{
  edm::Service<edm::RandomNumberGenerator> rng;
  rng_ = &rng->getEngine(lumi.index());
}

void cat::CATJetProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  runOnMC_ = !iEvent.isRealData();

  edm::Handle<pat::JetCollection> src;
  iEvent.getByToken(src_, src);

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  if ( !payloadName_.empty() ) {
    // temp measure - payloadName should be AK4PFchs, but PHYS14_25_V2 does not have uncertainty
    iSetup.get<JetCorrectionsRecord>().get(payloadName_,JetCorParColl);
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }

  JME::JetResolution jetResObj;
  JME::JetResolutionScaleFactor jetResSFObj;
  if ( runOnMC_ ) {
    jetResObj = JME::JetResolution(jetResFilePath_);
    jetResSFObj = JME::JetResolutionScaleFactor(jetResSFFilePath_);
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  const double rho = *rhoHandle;

  // for quark gluon likelihood calculation
  edm::Handle<edm::ValueMap<float>> qgHandle; 
  iEvent.getByToken(qgToken_, qgHandle);

  // for the different flavours
  std::vector<std::string> flavTagNames;
  std::vector<edm::Handle<edm::ValueMap<float>>> flavTagHandles;
  for ( auto token : flavTagTokens_ ) {
    flavTagHandles.push_back(edm::Handle<edm::ValueMap<float>>());
    iEvent.getByToken(token, flavTagHandles.back());
  }

  std::unique_ptr<cat::JetCollection>  out(new cat::JetCollection());

  for (auto aPatJetPointer = src->begin(); aPatJetPointer != src->end(); ++aPatJetPointer) {

    const pat::Jet& aPatJet = *aPatJetPointer;
    edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(src, aPatJetPointer - src->begin()));

    cat::Jet aJet(aPatJet);

    ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    float NHF = aPatJet.neutralHadronEnergyFraction();
    float NEMF = aPatJet.neutralEmEnergyFraction();
    float CHF = aPatJet.chargedHadronEnergyFraction();
    float MUF = aPatJet.muonEnergyFraction();
    float CEMF = aPatJet.chargedEmEnergyFraction();
    int NumConst = aPatJet.chargedMultiplicity()+aPatJet.neutralMultiplicity();
    int NumNeutralParticle =aPatJet.neutralMultiplicity();
    int CHM = aPatJet.chargedMultiplicity();
    float eta = aPatJet.eta();

    aJet.setChargedEmEnergyFraction(aPatJet.chargedEmEnergyFraction());
    bool looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7;
    bool tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7;
    bool tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) && abs(eta)<=2.7;

    if ( std::abs(eta) > 3.0 ) {
      looseJetID = (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 );
      tightJetID = (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 );
      tightLepVetoJetID = false;
    }
    else if (std::abs(eta) > 2.7){
      looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
      tightLepVetoJetID = false;
    }

    aJet.setLooseJetID( looseJetID );
    aJet.setTightJetID( tightJetID );
    aJet.setTightLepVetoJetID( tightLepVetoJetID );

    if( aPatJet.hasUserFloat("pileupJetId:fullDiscriminant") ) {
      aJet.setPileupJetId( aPatJet.userFloat("pileupJetId:fullDiscriminant") );
    }

    // aJet.addBDiscriminatorPair( aPatJet.bDiscriminator(btagNames_.at(0)) );

    if ( btagNames_.empty() ) {
      aJet.setBDiscriminators(aPatJet.getPairDiscri());
      // const std::vector<std::pair<std::string, float> > bpair = aPatJet.getPairDiscri();
      // for (unsigned int i =0; i < bpair.size(); i++){
      // 	cout << bpair[i].first <<endl;
      // }
    }
    else {
      for ( auto btagName : btagNames_ ) {
        aJet.addBDiscriminatorPair(std::make_pair(btagName, aPatJet.bDiscriminator(btagName) ));
      }
    }
    for ( int i=0, n=flavTagNames_.size(); i<n; ++i ) {
      const auto name = flavTagNames_[i];
      auto handle = flavTagHandles[i];
      const double value = handle.isValid() ? (*handle)[jetRef] : -999;
      aJet.addBDiscriminatorPair(std::make_pair(name, value));
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

    // calculate quark/gluon likelihood but only for AK4
    aJet.setQGLikelihood(-2.0);
    if ( qgHandle.isValid() ) {
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(src, aPatJetPointer - src->begin()));
      float qgLikelihood = (*qgHandle)[jetRef];
      aJet.setQGLikelihood(qgLikelihood);
      //aJet.setQGLikelihood(aPatJet.userFloat("QGTaggerAK4PFCHS:qgLikelihood"));
    }

    // setting JEC uncertainty
    if (!payloadName_.empty()){
      jecUnc->setJetEta(aJet.eta());
      jecUnc->setJetPt(aJet.pt()); // here you must use the CORRECTED jet pt
      double unc = jecUnc->getUncertainty(true);
      aJet.setShiftedEnUp( (1. + unc) );
      jecUnc->setJetEta(aJet.eta());
      jecUnc->setJetPt(aJet.pt()); // here you must use the CORRECTED jet pt
      unc = jecUnc->getUncertainty(false);
      aJet.setShiftedEnDown( (1. - unc) );
    }
    if (runOnMC_){
      // adding genJet
      auto genJet = aPatJet.genJetFwdRef();
      aJet.setGenJetRef(genJet);
      if (setGenParticle_) aJet.setGenParticleRef(aPatJet.genParticleRef());

      const double jetPt = aJet.pt();

      // Compute the JER
      JME::JetParameters jetPars = {{JME::Binning::JetPt, jetPt},
                                    {JME::Binning::JetEta, aJet.eta()},
                                    {JME::Binning::Rho, rho}};
      const double jetRes = jetResObj.getResolution(jetPars); // Note: this is relative resolution.
      const double cJER   = jetResSFObj.getScaleFactor(jetPars);
      const double cJERUp = jetResSFObj.getScaleFactor(jetPars, Variation::UP);
      const double cJERDn = jetResSFObj.getScaleFactor(jetPars, Variation::DOWN);

      // JER - apply scaling method if matched genJet is found,
      //       apply gaussian smearing method if unmatched
      if ( genJet.isNonnull() and deltaR(genJet->p4(), aJet.p4()) < 0.2
           and std::abs(genJet->pt()-jetPt) < jetRes*3*jetPt ) {
        const double genJetPt = genJet->pt();
        const double dPt = jetPt-genJetPt;
        const double fJER   = std::max(0., (genJetPt+dPt*cJER)/jetPt);
        const double fJERUp = std::max(0., (genJetPt+dPt*cJERUp)/jetPt);
        const double fJERDn = std::max(0., (genJetPt+dPt*cJERDn)/jetPt);

        aJet.setJER(fJER, fJERDn, fJERUp);
      }
      else {
        const double smear = CLHEP::RandGaussQ::shoot(rng_);

        const double fJER   = cJER   <= 1 ? 1 : 1+smear*jetRes*sqrt(cJER*cJER-1);
        const double fJERUp = cJERUp <= 1 ? 1 : 1+smear*jetRes*sqrt(cJERUp*cJERUp-1);
        const double fJERDn = cJERDn <= 1 ? 1 : 1+smear*jetRes*sqrt(cJERDn*cJERDn-1);

        aJet.setJER(fJER, fJERDn, fJERUp);
      }
    }

    out->push_back(aJet);
  }

  if (jecUnc) delete jecUnc;

  iEvent.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(CATJetProducer);
