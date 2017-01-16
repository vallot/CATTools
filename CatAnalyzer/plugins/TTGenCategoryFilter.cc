#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/DataFormats/interface/GenTop.h"

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
  enum InputType { IN_PartonTop, IN_PseudoTop, IN_GenTop } inputType_;

  bool vetoTau_;

  // For the parton top
  edm::EDGetTokenT<reco::GenParticleCollection> parton_srcToken_;
  edm::EDGetTokenT<int> parton_channelToken_;
  edm::EDGetTokenT<vint> parton_modesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> parton_jetToken_;
  const int nLepton_;

  // For the GenTop
  edm::EDGetTokenT<cat::GenTopCollection> genTop_srcToken_;
  int genTop_addJetCh_;

  // FIXME : classification by pseudotop, hadrons to be added

};

TTGenCategoryFilter::TTGenCategoryFilter(const edm::ParameterSet& pset):
  doInvert_(pset.getParameter<bool>("invert")),
  nLepton_(pset.getParameter<int>("nLepton"))
{
  const auto inputType = pset.getParameter<string>("inputType");
  if ( inputType == "PartonTop" ) inputType_ = IN_PartonTop;
  //else if ( inputType == "PseudoTop" ) inputType_ = IN_PseudoTop;
  else if ( inputType == "GenTop" ) inputType_ = IN_GenTop;
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
  else if ( inputType_ == IN_GenTop ) {
    vetoTau_ = pset.getParameter<bool>("vetoTau");
    genTop_srcToken_ = consumes<cat::GenTopCollection>(pset.getParameter<edm::InputTag>("src"));
    const auto addJetCh = pset.getParameter<std::string>("addJetChannel");
    genTop_addJetCh_ = 0;
    if      ( addJetCh == "TTBB" ) genTop_addJetCh_ = 1;
    else if ( addJetCh == "TTBJ" ) genTop_addJetCh_ = 2;
    else if ( addJetCh == "TTCC" ) genTop_addJetCh_ = 3;
    else if ( addJetCh == "TTJJ" ) genTop_addJetCh_ = 4;
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
    }; break;
    case IN_GenTop: {
      edm::Handle<cat::GenTopCollection> srcHandle;
      if ( !event.getByToken(genTop_srcToken_, srcHandle) ) break;
      if ( srcHandle->empty() ) break;

      const auto genTop = srcHandle->at(0);
      const int channelOption = vetoTau_ ? 0 : 1;
      if ( nLepton_ == 2 and !genTop.diLeptonic(channelOption) ) break;
      if ( nLepton_ == 1 and !genTop.semiLeptonic(channelOption) ) break;

      if ( genTop_addJetCh_ == 1 and genTop.NaddbJets20() <  2 ) break; // TTBB
      if ( genTop_addJetCh_ == 2 and genTop.NaddbJets20() != 1 ) break; // TTBJ
      if ( genTop_addJetCh_ == 3 and genTop.NaddcJets20() <  2 ) break; // TTCC
      if ( genTop_addJetCh_ == 4 and genTop.NaddJets20()  <  2 ) break; // TTJJ

      accept = true;
    }; break;
    case IN_PseudoTop:
      break;
  }

  if ( doInvert_ ) return !accept;
  return accept;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTGenCategoryFilter);
