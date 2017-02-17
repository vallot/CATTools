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
  else edm::LogError("TTGenCategoryFilter") << "Wrong input inputType. Choose among (PartonTop,GenTop)";

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
    genTop_addJetCh_ = 0; // Default value for "none of these category"
    if      ( addJetCh == "TTBB" ) genTop_addJetCh_ = 1;
    else if ( addJetCh == "TTBJ" ) genTop_addJetCh_ = 2;
    else if ( addJetCh == "TTCC" ) genTop_addJetCh_ = 3;
    else if ( addJetCh == "TTLF" ) genTop_addJetCh_ = 4;
  }

}

bool TTGenCategoryFilter::filter(edm::Event& event, const edm::EventSetup&)
{
  if ( inputType_ == IN_PartonTop ) {
    if ( nLepton_ < 0 ) return true;

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

    const bool accept = ( nLepton == nLepton_ );

    // return filter result, invert if "doInvert" is on
    //    if ( doInvert_ ) return !accept;
    //    return accept;
    // .. in shorter syntax:
    //   doInvert accpet -> result
    //      T       T    ->   F
    //      T       F    ->   T
    //      F       T    ->   T
    //      F       F    ->   F
    return doInvert_ xor accept;
  }
  else if ( inputType_ == IN_GenTop ) {
    edm::Handle<cat::GenTopCollection> srcHandle;
    event.getByToken(genTop_srcToken_, srcHandle);
    if ( srcHandle->empty() ) return false;

    const auto genTop = srcHandle->at(0);
    const int channelOption = vetoTau_ ? 0 : 1;

    // Same logic with the above
    const bool acceptCh = ( nLepton_ == 2 and genTop.diLeptonic(channelOption) ) or
                          ( nLepton_ == 1 and genTop.semiLeptonic(channelOption) );
    if ( nLepton_ >= 0 and (doInvert_ xor acceptCh) == false ) return false;

    // non-ttjj cases
    if ( genTop.NaddJets20() < 2 ) return (genTop_addJetCh_ == 0);

    if      ( genTop.NaddbJets20() >= 2 and genTop_addJetCh_ == 1 ) return true; // TTBB (2 b jets)
    else if ( genTop.NaddbJets20() == 1 and genTop_addJetCh_ == 2 ) return true; // TTBJ (1 b jet)
    else if ( genTop.NaddcJets20() >= 2 and genTop_addJetCh_ == 3 ) return true; // TTCC (2 c jets)
    return genTop_addJetCh_ == 4; // others are TTLF
  }

  return false; // dummy return
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTGenCategoryFilter);
