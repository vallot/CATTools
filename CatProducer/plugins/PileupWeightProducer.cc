#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/View.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class PileupWeightProducer : public edm::stream::EDProducer<>
{
public:
  PileupWeightProducer(const edm::ParameterSet& pset);
  ~PileupWeightProducer() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::LumiReWeighting lumiWeights_, lumiWeightsUp_, lumiWeightsDn_;

  const bool isStandardWeight_;
  typedef std::vector<PileupSummaryInfo> PUInfos;
  edm::EDGetTokenT<PUInfos> puToken_;

  std::vector<double> simpleWeights_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

};

PileupWeightProducer::PileupWeightProducer(const edm::ParameterSet& pset):
  isStandardWeight_(pset.getParameter<bool>("isStandardWeight")),
  puToken_(consumes<PUInfos>(pset.getParameter<edm::InputTag>("pileupInfo")))
{
  if ( isStandardWeight_ )
  {
    std::vector<double> pileupMC = pset.getParameter<std::vector<double> >("pileupMC");
    std::vector<double> pileupRD = pset.getParameter<std::vector<double> >("pileupRD");
    std::vector<double> pileupUp = pset.getParameter<std::vector<double> >("pileupUp");
    std::vector<double> pileupDn = pset.getParameter<std::vector<double> >("pileupDn");
    const double sumWMC = std::accumulate(pileupMC.begin(), pileupMC.end(), 0.);
    const double sumWRD = std::accumulate(pileupRD.begin(), pileupRD.end(), 0.);
    const double sumWUp = std::accumulate(pileupUp.begin(), pileupUp.end(), 0.);
    const double sumWDn = std::accumulate(pileupDn.begin(), pileupDn.end(), 0.);

    std::vector<float> pileupMCTmp;
    std::vector<float> pileupRDTmp;
    std::vector<float> pileupUpTmp, pileupDnTmp;
    for ( int i=0, n=min(pileupMC.size(), pileupRD.size()); i<n; ++i )
    {
      pileupMCTmp.push_back(pileupMC[i]/sumWMC);
      pileupRDTmp.push_back(pileupRD[i]/sumWRD);
      pileupUpTmp.push_back(pileupUp[i]/sumWUp);
      pileupDnTmp.push_back(pileupDn[i]/sumWDn);
    }
    lumiWeights_ = edm::LumiReWeighting(pileupMCTmp, pileupRDTmp);
    lumiWeightsUp_ = edm::LumiReWeighting(pileupMCTmp, pileupUpTmp);
    lumiWeightsDn_ = edm::LumiReWeighting(pileupMCTmp, pileupDnTmp);
  }
  else
  {
    std::cerr << "!!PileupWeightProducer!! We are using NON STANDARD method for the pileup reweight.\n"
              << "                         This weight values are directly from reco vertex\n";
    vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));
    simpleWeights_ = pset.getParameter<std::vector<double> >("simpleWeights");
    const double sumW = std::accumulate(simpleWeights_.begin(), simpleWeights_.end(), 0.);
    if ( sumW > 0 ) { for ( auto& w : simpleWeights_ ) { w /= sumW; } }
  }

  produces<int>("nTrueInteraction");
  produces<float>("");
  produces<float>("up");
  produces<float>("dn");

}

void PileupWeightProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<int> nTrueIntr(new int(-1));
  std::auto_ptr<float> weight(new float(1.));
  std::auto_ptr<float> weightUp(new float(1.));
  std::auto_ptr<float> weightDn(new float(1.));

  if ( !event.isRealData() ){
    if ( isStandardWeight_ ){
      edm::Handle<std::vector<PileupSummaryInfo> > puHandle;
      event.getByToken(puToken_, puHandle);

      for ( auto& puInfo : *puHandle ){
        //const int nIntr = puInfo.getPU_NumInteractions();
        const int bx = puInfo.getBunchCrossing();

        if ( bx == 0 ){
          *nTrueIntr = puInfo.getTrueNumInteractions();
          break;
        }
      }

      if ( *nTrueIntr > 0 ) {
        *weight   = lumiWeights_.weight(*nTrueIntr);
        *weightUp = lumiWeightsUp_.weight(*nTrueIntr);
        *weightDn = lumiWeightsDn_.weight(*nTrueIntr);
      }
    }
    else {
      edm::Handle<reco::VertexCollection> vertexHandle;
      event.getByToken(vertexToken_, vertexHandle);

      const int nPVBin = std::min(simpleWeights_.size(), vertexHandle->size()) - 1;
      if ( nPVBin >= 0 ){
        *weight   = simpleWeights_[nPVBin];
        *weightUp = simpleWeights_[nPVBin];
        *weightDn = simpleWeights_[nPVBin];
      }
    }
  }

  event.put(nTrueIntr, "nTrueInteraction");
  event.put(weight  , "");
  event.put(weightUp, "up");
  event.put(weightDn, "dn");
}

DEFINE_FWK_MODULE(PileupWeightProducer);
