#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/Trigger.h"

#include <memory>
#include <vector>
#include <string>

class CATSkimEventFilter : public edm::stream::EDFilter<>
{
public:
  CATSkimEventFilter(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<cat::ElectronCollection> electronsToken_;
  edm::EDGetTokenT<cat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<cat::JetCollection> jetsToken_;
  
  const double minLeptonPt_, maxLeptonAbseta_;
  const unsigned int minNLeptons_;
  const std::vector<std::string> electronIdNames_;

  const double minJetPt_, maxJetAbseta_;
  const unsigned int minNJets_;
};

CATSkimEventFilter::CATSkimEventFilter(const edm::ParameterSet& pset):
  minLeptonPt_(pset.getParameter<double>("minLeptonPt")),
  maxLeptonAbseta_(pset.getParameter<double>("maxLeptonAbseta")),
  minNLeptons_(pset.getParameter<unsigned int>("minNLeptons")),
  electronIdNames_(pset.getParameter<std::vector<std::string>>("electronIdNames")),
  minJetPt_(pset.getParameter<double>("minJetPt")),
  maxJetAbseta_(pset.getParameter<double>("maxJetAbseta")),
  minNJets_(pset.getParameter<unsigned int>("minNJets"))
{
  electronsToken_ = consumes<cat::ElectronCollection>(pset.getParameter<edm::InputTag>("electrons"));
  muonsToken_ = consumes<cat::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));
  jetsToken_ = consumes<cat::JetCollection>(pset.getParameter<edm::InputTag>("jets"));
}

using namespace std;

bool CATSkimEventFilter::filter(edm::Event& event, const edm::EventSetup&)
{
  if ( minNLeptons_ == 0 and minNJets_ == 0 ) return true;

  edm::Handle<cat::ElectronCollection> electronsHandle;
  event.getByToken(electronsToken_, electronsHandle);

  edm::Handle<cat::MuonCollection> muonsHandle;
  event.getByToken(muonsToken_, muonsHandle);

  edm::Handle<cat::JetCollection> jetsHandle;
  event.getByToken(jetsToken_, jetsHandle);

  unsigned int nLeptons = 0, nJets = 0;

  for ( auto& ele : *electronsHandle ) {
    const double scAbseta = std::abs(ele.scEta());
    const double abseta = std::abs(ele.eta());
    if ( std::min(scAbseta, abseta) > maxLeptonAbseta_ ) continue;

    const double pt = ele.pt()*ele.shiftedEnDown()*std::min(1.F, ele.smearedScale());
    if ( pt < minLeptonPt_ ) continue;

    bool passId = false; //(ele.isVeto() or !ele.isMediumMVA() or !ele.isTightMVA());
    passId = (passId or ele.snuID());
    for ( auto& idName : electronIdNames_ ) {
      if ( ele.electronID(idName) != 0 ) {
        passId = true;
        break;
      }
    }
    if ( !passId ) continue;

    ++nLeptons;
  }

  for ( auto& mu : *muonsHandle ) {
    const double abseta = std::abs(mu.eta());
    if ( abseta > maxLeptonAbseta_ ) continue;

    const double pt = mu.pt()*mu.shiftedEnDown();
    if ( pt < minLeptonPt_ ) continue;

    bool passId = (mu.isLooseMuon() or mu.isMediumMuon() or mu.isTightMuon());

    if ( !passId ) continue;
    ++nLeptons;
  }

  for ( auto& jet : *jetsHandle ) {
    const double abseta = std::abs(jet.eta());
    if ( abseta > maxLeptonAbseta_ ) continue;

    const double pt = jet.pt()*jet.shiftedEnDown()*jet.smearedResDown();
    if ( pt < minJetPt_ ) continue;

    bool passId = (jet.looseJetID() or jet.tightJetID() or jet.tightLepVetoJetID() );
    
    if ( !passId ) continue;
    ++nJets;
  }

  if ( nLeptons < minNLeptons_ ) return false;
  if ( nJets < minNJets_ ) return false;

  return true;
}

DEFINE_FWK_MODULE(CATSkimEventFilter);

