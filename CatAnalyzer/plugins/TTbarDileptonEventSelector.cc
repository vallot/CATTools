#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "CommonTools/Utils/interface/PtComparator.h"

using namespace std;

namespace cat {

class TTbarDileptonEventSelector : public edm::stream::EDProducer<>
{
public:
  TTbarDileptonEventSelector(const edm::ParameterSet& pset);
  void produce(edm::Event & event, const edm::EventSetup&) override;

private:
  typedef cat::Muon TMuon;
  typedef cat::Electron TElectron;;
  typedef cat::Jet TJet;
  typedef cat::MET TMET;
  edm::EDGetTokenT<edm::View<TMuon> > muonToken_;
  edm::EDGetTokenT<edm::View<TElectron> > electronToken_;
  edm::EDGetTokenT<edm::View<TJet> > jetToken_;
  edm::EDGetTokenT<edm::View<TMET> > metToken_;

private:
  bool isGoodMuon(const TMuon& mu)
  {
    if ( mu.pt() <= 20 or std::abs(mu.eta()) >= 2.4 ) return false;
    if ( mu.relIso(0.4) >= 0.12 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isVetoMuon(const TMuon& mu)
  {
    if ( mu.pt() <= 20 or std::abs(mu.eta()) >= 2.4 ) return false;
    if ( mu.relIso(0.4) >= 0.2 ) return false;
    if ( !mu.isLooseMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const TElectron& el)
  {
    if ( el.pt() <= 20 or std::abs(el.eta()) >= 2.4 ) return false;
    //if ( el.relIso(0.3) >= 0.11 ) return false;
    if ( !el.electronID(eleIdName_) ) return false;
    if ( !el.isPF() or !el.passConversionVeto() ) return false;
    const double scEta = std::abs(el.scEta());
    if ( scEta >= 1.4442 and scEta <= 1.566 ) return false;
    return true;
  }
  bool isVetoElectron(const TElectron& el)
  {
    if ( el.pt()  <= 20 or std::abs(el.eta()) >= 2.4 ) return false;
    if ( !el.electronID(eleVetoIdName_) ) return false;
    const double scEta = std::abs(el.scEta());
    if ( scEta >= 1.4442 and scEta <= 1.566 ) return false;
    const double relIso03 = el.relIso(0.3);
    if ( scEta <= 1.4442 and relIso03 >= 0.1649 ) return false;
    if ( scEta >= 1.566  and relIso03 >= 0.2075 ) return false;
    return true;
  }
  bool isBjet(const TJet& jet)
  {
    if ( jet.bDiscriminator(bTagName_) >= 0.814 ) return true;
    return false;
  }

private:
  typedef reco::Candidate::LorentzVector LorentzVector;
  typedef reco::CompositePtrCandidate CRCand;
  typedef std::vector<CRCand> CRCandColl;
  typedef std::vector<double> doubles;
  enum CHANNEL {
    CH_NONE=0, CH_MUMU, CH_ELEL, CH_MUEL
  };

  const bool keepVetoLeptons_, checkOverlapFromVetoLeptons_;
  const bool sortByBtag_;
  const std::string eleIdName_, eleVetoIdName_;
  const std::string bTagName_;

};

}

using namespace cat;

TTbarDileptonEventSelector::TTbarDileptonEventSelector(const edm::ParameterSet& pset):
  keepVetoLeptons_(pset.getParameter<bool>("keepVetoLeptons")),
  checkOverlapFromVetoLeptons_(pset.getParameter<bool>("checkOverlapFromVetoLepton")),
  sortByBtag_(pset.getParameter<bool>("sortByBtag")),
  eleIdName_(pset.getParameter<std::string>("eleIdName")),
  eleVetoIdName_(pset.getParameter<std::string>("eleVetoIdName")),
  bTagName_(pset.getParameter<std::string>("bTagName"))
{
  muonToken_ = consumes<edm::View<TMuon> >(pset.getParameter<edm::InputTag>("muons"));
  electronToken_ = consumes<edm::View<TElectron> >(pset.getParameter<edm::InputTag>("electrons"));
  jetToken_ = consumes<edm::View<TJet> >(pset.getParameter<edm::InputTag>("jets"));
  metToken_ = consumes<edm::View<TMET> >(pset.getParameter<edm::InputTag>("mets"));

  produces<int>("channel");
  produces<int>("nLeptons");
  produces<int>("nVetoLeptons");
  produces<int>("nBjets");
  produces<float>("met");
  produces<float>("metphi");
  produces<std::vector<reco::CandidatePtr> >("leptons");
  produces<std::vector<reco::CandidatePtr> >("jets");
}

void TTbarDileptonEventSelector::produce(edm::Event& event, const edm::EventSetup&)
{
  // Output contents
  int channel = CH_NONE, nLeptons = 0, nVetoLeptons = 0;
  int nBjets = 0;
  std::auto_ptr<std::vector<reco::CandidatePtr> > out_leptons(new std::vector<reco::CandidatePtr>);
  std::auto_ptr<std::vector<reco::CandidatePtr> > out_jets(new std::vector<reco::CandidatePtr>);

  // MET first
  edm::Handle<edm::View<TMET> > metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& metP4 = metHandle->at(0).p4();

  do
  {
    // Do the leptons
    edm::Handle<edm::View<TMuon> > muonHandle;
    event.getByToken(muonToken_, muonHandle);

    edm::Handle<edm::View<TElectron> > electronHandle;
    event.getByToken(electronToken_, electronHandle);

    std::vector<reco::CandidatePtr> leptons;
    std::set<reco::CandidatePtr> vetoLeptons;
    for ( int i=0, n=muonHandle->size(); i<n; ++i )
    {
      auto& p = muonHandle->at(i);
      if ( p.pt() < 20 or abs(p.eta()) >= 2.4 ) continue;

      reco::CandidatePtr muonPtr = reco::CandidatePtr(muonHandle, i);

      if ( isGoodMuon(p) ) leptons.push_back(muonPtr);
      if ( isVetoMuon(p) ) vetoLeptons.insert(muonPtr);
    }
    for ( int i=0, n=electronHandle->size(); i<n; ++i )
    {
      auto& p = electronHandle->at(i);
      if ( p.pt() < 20 or std::abs(p.eta()) > 2.4 ) continue;

      reco::CandidatePtr electronPtr = reco::CandidatePtr(electronHandle, i);

      if ( isGoodElectron(p) ) leptons.push_back(electronPtr);
      if ( isVetoElectron(p) ) vetoLeptons.insert(electronPtr);
    }

    // Partial sort to select leading 2 leptons
    auto GtByPtPtr = [](reco::CandidatePtr a, reco::CandidatePtr b){return a->pt() > b->pt();};
    std::nth_element(leptons.begin(), leptons.begin()+2, leptons.end(), GtByPtPtr);

    nLeptons = leptons.size();
    // Remove selected leptons from veto lepton collection
    if ( nLeptons > 0 ) vetoLeptons.erase(leptons.at(0));
    if ( nLeptons > 1 ) vetoLeptons.erase(leptons.at(1));
    nVetoLeptons = vetoLeptons.size();

    // At least two selected leptons
    if ( nLeptons < 2 ) break;
    auto lepton1 = leptons.at(0), lepton2 = leptons.at(1);

    // Set channel
    const int pdgId1 = std::abs(lepton1->pdgId());
    const int pdgId2 = std::abs(lepton2->pdgId());
    if ( pdgId1 == 13 and pdgId2 == 13 ) channel = CH_MUMU;
    else if ( pdgId1 == 11 and pdgId2 == 11 ) channel = CH_ELEL;
    else channel = CH_MUEL;

     // Minimal dilepton mass cut
    if ( (lepton1->p4()+lepton2->p4()).mass() < 20 ) break;

    // Put selected leptons into the output collection
    out_leptons->push_back(lepton1);
    out_leptons->push_back(lepton2);

    if ( keepVetoLeptons_ )
    {
      // Append veto leptons to the output lepton collection
      for ( auto x : vetoLeptons ) out_leptons->push_back(x);
      // Sort veto leptons keeping selected leptons front
      std::sort(out_leptons->begin()+2, out_leptons->end(),GtByPtPtr);
    }

    // Continue to jets
    edm::Handle<edm::View<TJet> > jetHandle;
    event.getByToken(jetToken_, jetHandle);
    for ( int i=0, n=jetHandle->size(); i<n; ++i )
    {
      auto& jet = jetHandle->at(i);

      if ( jet.pt() < 30 or std::abs(jet.eta()) > 2.5 ) continue;
      if ( deltaR(jet.p4(), lepton1->p4()) < 0.5 ) continue;
      if ( deltaR(jet.p4(), lepton2->p4()) < 0.5 ) continue;
      if ( checkOverlapFromVetoLeptons_ )
      {
        bool isOverlap = false;
        for ( auto& l : vetoLeptons ) if ( deltaR(jet.p4(), l->p4()) < 0.5 ) { isOverlap = true; break; }
        if ( isOverlap ) continue;
      }

      out_jets->push_back(reco::CandidatePtr(jetHandle, i));
      if ( isBjet(jet) ) ++nBjets;
    }

    // Sort by b-discriminator (up to nBjet, to be used in ttbb analysis)
    if ( sortByBtag_ )
    {
      const string bTagName(bTagName_);
      std::nth_element(out_jets->begin(), out_jets->begin()+nBjets, out_jets->end(),
                       [bTagName](reco::CandidatePtr a, reco::CandidatePtr b){
                         const double bTag1 = dynamic_cast<const TJet&>(*a).bDiscriminator(bTagName);
                         const double bTag2 = dynamic_cast<const TJet&>(*b).bDiscriminator(bTagName);
                         return bTag1 > bTag2;});
      // Sort by pt for the others
      std::sort(out_jets->begin()+nBjets, out_jets->end(), GtByPtPtr);
    }
    else
    {
      std::sort(out_jets->begin(), out_jets->end(), GtByPtPtr);
    }
  } while ( false );

  event.put(std::auto_ptr<int>(new int(channel)), "channel");
  event.put(std::auto_ptr<int>(new int(nLeptons)), "nLeptons");
  event.put(std::auto_ptr<int>(new int(nVetoLeptons)), "nVetoLeptons");
  event.put(std::auto_ptr<int>(new int(nBjets)), "nBjets");
  event.put(std::auto_ptr<float>(new float(metP4.pt())), "met");
  event.put(std::auto_ptr<float>(new float(metP4.phi())), "metphi");
  event.put(out_leptons, "leptons");
  event.put(out_jets, "jets");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTbarDileptonEventSelector);

