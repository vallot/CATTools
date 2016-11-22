#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"

using namespace std;
using namespace cat;

class TTGenCategoryFilter : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TTGenCategoryFilter(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;

private:
  const bool doInvert_;

  typedef std::vector<int> vint;
  enum InputType { IN_PartonTop, IN_PseudoTop, IN_Hadron } inputType_;

  // For the parton top
  edm::EDGetTokenT<reco::GenParticleCollection> parton_srcToken_;
  edm::EDGetTokenT<int> parton_channelToken_;
  edm::EDGetTokenT<vint> parton_modesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> parton_jetToken_;
  const int nLepton_;
  bool vetoTau_;

  // FIXME : classification by pseudotop, hadrons to be added

};

TTGenCategoryFilter::TTGenCategoryFilter(const edm::ParameterSet& pset):
  doInvert_(pset.getParameter<bool>("invert")),
  nLepton_(pset.getParameter<int>("nLepton"))
{
  const auto inputType = pset.getParameter<string>("inputType");
  if ( inputType == "PartonTop" ) inputType_ = IN_PartonTop;
  //else if ( inputType == "PseudoTop" ) inputType_ = IN_PseudoTop;
  //else if ( inputType == "Hadron" ) inputType_ = IN_Hadron;
  else edm::LogError("TTGenCategoryFilter") << "Wrong input inputType. Choose among (PartonTop,)";

  const auto addJetType = pset.getParameter<string>("addJetType");

  if ( inputType_ == IN_PartonTop ) {
    vetoTau_ = pset.getParameter<bool>("vetoTau");
    const auto label = pset.getParameter<edm::InputTag>("src");
    const auto labelName = label.label();
    parton_srcToken_ = consumes<reco::GenParticleCollection>(label);
    parton_channelToken_ = consumes<int>(edm::InputTag(labelName, "channel"));
    parton_modesToken_ = consumes<vint>(edm::InputTag(labelName, "modes"));
    parton_jetToken_ = consumes<reco::GenJetCollection>(edm::InputTag(labelName, "qcdJets"));
  }

}

bool TTGenCategoryFilter::filter(edm::Event& event, const edm::EventSetup&)
{
  bool accept = false;
  switch ( inputType_ ) {
    case IN_PartonTop: {
      edm::Handle<reco::GenParticleCollection> srcHandle;
      edm::Handle<int> channelHandle;
      edm::Handle<vint> modesHandle;
      edm::Handle<reco::GenJetCollection> jetHandle;

      event.getByToken(parton_srcToken_, srcHandle);
      event.getByToken(parton_channelToken_, channelHandle);
      event.getByToken(parton_modesToken_, modesHandle);
      event.getByToken(parton_jetToken_, jetHandle);

      //const int channel = *channelHandle;
      const int mode1 = modesHandle->at(0);
      const int mode2 = modesHandle->at(1);

      int nLepton = 0;
      if ( mode1%3 != 0 and (!vetoTau_ or mode1 < CH_TAU_HADRON) ) ++nLepton;
      if ( mode2%3 != 0 and (!vetoTau_ or mode2 < CH_TAU_HADRON) ) ++nLepton;

      if ( nLepton != nLepton_ ) break;

      accept = true;
    }
    case IN_PseudoTop:
    case IN_Hadron:
      break;
  }

  if ( doInvert_ ) return !accept;
  return accept;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTGenCategoryFilter);
