#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <vector>
#include <string>

class EventWeightThresholdFilter : public edm::stream::EDFilter<>
{
public:
  EventWeightThresholdFilter(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;

private:
  //typedef std::vector<bool> vbool;
  //typedef std::vector<int> vint;
  typedef std::vector<std::string> strings;
  bool combineByOr_;

  const bool isLessThan_;
  const double threshold_;
  edm::EDGetTokenT<double> weightToken_;
  edm::EDGetTokenT<float> weightTokenF_;
  edm::EDGetTokenT<int> weightTokenI_;

  enum class TYPE { D, F, I } type_;

};

EventWeightThresholdFilter::EventWeightThresholdFilter(const edm::ParameterSet& pset):
  isLessThan_(pset.getParameter<bool>("isLessThan")),
  threshold_(pset.getParameter<double>("threshold"))
{
  const auto typeStr = pset.getParameter<std::string>("type");
  if ( typeStr == "double" ) {
    type_ = TYPE::D;
    weightToken_ = consumes<double>(pset.getParameter<edm::InputTag>("src"));
  }
  else if ( typeStr == "float" ) {
    type_ = TYPE::F;
    weightTokenF_ = consumes<float>(pset.getParameter<edm::InputTag>("src"));
  }
  else if ( typeStr == "int" ) {
    type_ = TYPE::I;
    weightTokenI_ = consumes<int>(pset.getParameter<edm::InputTag>("src"));
  }
  else edm::LogError("EventWeightThresholdFilter") << "Wrong input to \"type\" option, it was \"" << typeStr << ".\n"
                                              << "This should be chosen among (\"double\", \"float\", \"int\")\n";

}

bool EventWeightThresholdFilter::filter(edm::Event& event, const edm::EventSetup&)
{
  using namespace std;

  double w = -1;
  if ( type_ == TYPE::D ) {
    edm::Handle<double> handle;
    event.getByToken(weightToken_, handle);
    w = *handle;
  }
  else if ( type_ == TYPE::F ) {
    edm::Handle<float> handle;
    event.getByToken(weightTokenF_, handle);
    w = *handle;
  }
  else if ( type_ == TYPE::I ) {
    edm::Handle<int> handle;
    event.getByToken(weightTokenI_, handle);
    w = *handle;
  }

  if ( isLessThan_ ) return w < threshold_;
  return w >= threshold_;
}

DEFINE_FWK_MODULE(EventWeightThresholdFilter);
