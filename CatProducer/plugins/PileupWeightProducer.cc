#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Common/interface/View.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class PileupWeightProducer : public edm::EDProducer
{
public:
  PileupWeightProducer(const edm::ParameterSet& pset);
  ~PileupWeightProducer() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::LumiReWeighting lumiWeights_, lumiWeightsUp_, lumiWeightsDn_;

  typedef std::vector<PileupSummaryInfo> PUInfos;
  edm::EDGetTokenT<PUInfos> puToken_;

};

PileupWeightProducer::PileupWeightProducer(const edm::ParameterSet& pset)
{
  std::vector<double> pileupMC = pset.getParameter<std::vector<double> >("pileupMC");
  std::vector<double> pileupRD = pset.getParameter<std::vector<double> >("pileupRD");
  std::vector<double> pileupUp = pset.getParameter<std::vector<double> >("pileupUp");
  std::vector<double> pileupDn = pset.getParameter<std::vector<double> >("pileupDn");

  std::vector<float> pileupMCTmp;
  std::vector<float> pileupRDTmp;
  std::vector<float> pileupUpTmp, pileupDnTmp;
  for ( int i=0, n=min(pileupMC.size(), pileupRD.size()); i<n; ++i )
  {
    pileupMCTmp.push_back(pileupMC[i]);
    pileupRDTmp.push_back(pileupRD[i]);
    pileupUpTmp.push_back(pileupUp[i]);
    pileupDnTmp.push_back(pileupDn[i]);
  }
  lumiWeights_ = edm::LumiReWeighting(pileupMCTmp, pileupRDTmp);
  lumiWeightsUp_ = edm::LumiReWeighting(pileupMCTmp, pileupUpTmp);
  lumiWeightsDn_ = edm::LumiReWeighting(pileupMCTmp, pileupDnTmp);

  puToken_ = consumes<PUInfos>(edm::InputTag("addPileupInfo"));

  produces<int>("nTrueInteraction");
  produces<double>("");
  produces<double>("up");
  produces<double>("dn");

}

void PileupWeightProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<std::vector<PileupSummaryInfo> > puHandle;
  event.getByToken(puToken_, puHandle);

  std::auto_ptr<int> nTrueIntr(new int(-1));
  std::auto_ptr<double> weight(new double(1.));
  std::auto_ptr<double> weightUp(new double(1.));
  std::auto_ptr<double> weightDn(new double(1.));

  if ( !event.isRealData() and puHandle.isValid() )
  {
    for ( auto& puInfo : *puHandle )
    {
      //const int nIntr = puInfo.getPU_NumInteractions();
      const int bx = puInfo.getBunchCrossing();

      if ( bx == 0 )
      {
        *nTrueIntr = puInfo.getTrueNumInteractions();
      }
    }
    
    if ( *nTrueIntr > 0 )
    {
      *weight   = lumiWeights_.weight(*nTrueIntr);
      *weightUp = lumiWeightsUp_.weight(*nTrueIntr);
      *weightDn = lumiWeightsDn_.weight(*nTrueIntr);
    }
  }

  event.put(nTrueIntr, "nTrueInteraction");
  event.put(weight  , "");
  event.put(weightUp, "up");
  event.put(weightDn, "dn");
}

DEFINE_FWK_MODULE(PileupWeightProducer);
