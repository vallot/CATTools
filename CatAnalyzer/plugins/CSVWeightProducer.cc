#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/CatAnalyzer/interface/CSVHelper.h"

#include <vector>
#include <memory>

class CSVWeightProducer : public edm::stream::EDProducer<>
{
public:
  CSVWeightProducer(const edm::ParameterSet& pset);
  ~CSVWeightProducer() {}
  void produce(edm::Event& event, const edm::EventSetup&) override;

private:
  const double minPt_, maxEta_;
  typedef std::vector<float> vfloat;
  edm::EDGetTokenT<edm::View<reco::Candidate>> lepToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;

  std::unique_ptr<CSVHelper> helper_;
};

CSVWeightProducer::CSVWeightProducer(const edm::ParameterSet& pset):
  minPt_(pset.getParameter<double>("minPt")),
  maxEta_(pset.getParameter<double>("maxEta"))
{
  jetToken_ = consumes<cat::JetCollection>(pset.getParameter<edm::InputTag>("src"));
  lepToken_ = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("leptons"));
  auto csvFileNameHF = pset.getParameter<std::string>("csvSFHF");
  auto csvFileNameLF = pset.getParameter<std::string>("csvSFLF");
  auto csvFileHF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileNameHF).fullPath();
  auto csvFileLF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileNameLF).fullPath();

  helper_.reset(new CSVHelper(csvFileHF, csvFileLF));

  produces<float>();
  produces<vfloat>("syst");
}

void CSVWeightProducer::produce(edm::Event& event, const edm::EventSetup&)
{
  using namespace std;

  std::unique_ptr<float> weight(new float(1.));
  std::unique_ptr<vfloat> weightSysts(new vfloat());

  if ( !event.isRealData() ) {
    edm::Handle<cat::JetCollection> jetHandle;
    event.getByToken(jetToken_, jetHandle);

    edm::Handle<edm::View<reco::Candidate>> lepHandle;
    event.getByToken(lepToken_, lepHandle);

    cat::JetCollection jets;
    for ( auto& jet : *jetHandle) {
      if ( jet.pt() < minPt_ or std::abs(jet.eta()) > maxEta_ ) continue;
      const bool isOverlap = [&](){
        if ( !lepHandle.isValid() ) return false;
        for ( auto& lep : *lepHandle ) {
          if ( deltaR(lep, jet) < 0.3 ) return true;
        }
        return false;
      }();
      if ( isOverlap ) continue;

      jets.push_back(jet);
    }

    *weight = helper_->getCSVWeight(jets, 0);
    for ( int i=7; i<25; ++i ) {
      weightSysts->push_back(helper_->getCSVWeight(jets, i));
    }
  }

  event.put(std::move(weight));
  event.put(std::move(weightSysts), "syst");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CSVWeightProducer);
