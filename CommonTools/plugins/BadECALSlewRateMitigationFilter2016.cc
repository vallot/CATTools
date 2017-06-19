#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/EDCollection.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class BadECALSlewRateMitigationFilter2016 : public edm::stream::EDFilter<>
{
public:
  BadECALSlewRateMitigationFilter2016(const edm::ParameterSet& pset);

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<bool> dupECALClusterToken_;
  edm::EDGetTokenT<edm::EDCollection<DetId>> hitsNotReplacedToken_;

  const bool doFilter_;
};

BadECALSlewRateMitigationFilter2016::BadECALSlewRateMitigationFilter2016(const edm::ParameterSet& pset):
  doFilter_(pset.getUntrackedParameter<bool>("doFilter", false))
{
  dupECALClusterToken_ = consumes<bool>(edm::InputTag("particleFlowEGammaGSFixed:dupECALClusters"));
  hitsNotReplacedToken_ = consumes<edm::EDCollection<DetId>>(edm::InputTag("ecalMultiAndGSGlobalRecHitEB:hitsNotReplaced"));
  
  produces<bool>();
}

bool BadECALSlewRateMitigationFilter2016::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  bool isAccepted = true;

  do {
    if ( !event.isRealData() ) break; // always accept MC

    edm::Handle<bool> dupECALClusterHandle;
    event.getByToken(dupECALClusterToken_, dupECALClusterHandle);

    edm::Handle<edm::EDCollection<DetId>> hitsNotReplacedHandle;
    event.getByToken(hitsNotReplacedToken_, hitsNotReplacedHandle);

    if ( dupECALClusterHandle.isValid() and *dupECALClusterHandle == true ) {
      // Duplicated ECAL clusters must be false
      isAccepted = false;
      break;
    }

    if ( hitsNotReplacedHandle.isValid() and !hitsNotReplacedHandle->empty() ) {
      // There should no "hits not replaced"
      isAccepted = false;
      break;
    }

  } while (false);

  event.put(std::move(std::make_unique<bool>(isAccepted)));

  if ( !doFilter_ ) return true;
  return isAccepted;
}

DEFINE_FWK_MODULE(BadECALSlewRateMitigationFilter2016);
