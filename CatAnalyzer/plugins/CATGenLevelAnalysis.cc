#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"

#include <iostream>

using namespace std;

class CATGenLevelAnalysis : public edm::EDAnalyzer
{
public:
  CATGenLevelAnalysis(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<int> channelToken_;
  edm::EDGetTokenT<std::vector<int> > modesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonsToken_;

  TH1F* hLepton1Pt_;
  TH1F* hLepton1Eta_;
  TH1F* hLepton1Phi_;
  TH1F* hLepton2Pt_;
  TH1F* hLepton2Eta_;
  TH1F* hLepton2Phi_;

  TH1F* hJet1Pt_;
  TH1F* hJet1Eta_;
  TH1F* hJet1Phi_;

  TH1F* hLeptonJetDeltaR_;

};

//using namespace cat;

CATGenLevelAnalysis::CATGenLevelAnalysis(const edm::ParameterSet& pset)
{
  partonsToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("partons"));
  channelToken_ = consumes<int>(pset.getParameter<edm::InputTag>("channel"));
  modesToken_ = consumes<std::vector<int> >(pset.getParameter<edm::InputTag>("modes"));

  edm::Service<TFileService> fs;
  hLepton1Pt_ = fs->make<TH1F>("hLepton1Pt", "lepton 1 p_{T};p_{T} (GeV/c);Events", 100, 0, 100);
  // Do the same for the other histograms
}

void CATGenLevelAnalysis::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<int> channelHandle;
  event.getByToken(channelToken_, channelHandle);

  edm::Handle<std::vector<int> > modesHandle;
  event.getByToken(modesToken_, modesHandle);

  edm::Handle<reco::GenParticleCollection> partonsHandle;
  event.getByToken(partonsToken_, partonsHandle);

  cout << "CHANNEL = " << *channelHandle << endl;
  for ( int i=0, n=modesHandle->size(); i<n; ++i )
  {
    cout << "MODE" << i << " = " << modesHandle->at(i) << endl;
  }

  for ( int i=0, n=partonsHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = partonsHandle->at(i);

    cout << p.pt() << ' ' << p.pdgId() << ' ' << p.pt() << ' ' << p.eta() << ' ' << p.phi() << endl;
  }
}

DEFINE_FWK_MODULE(CATGenLevelAnalysis);

