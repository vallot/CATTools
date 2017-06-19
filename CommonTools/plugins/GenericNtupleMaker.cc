#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "CATTools/CommonTools/interface/Consumers.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;
using namespace edm;

using namespace cat;

class GenericNtupleMaker : public edm::EDAnalyzer
{
public:
  GenericNtupleMaker(const edm::ParameterSet& pset);

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup) override;

private:
  typedef edm::ParameterSet PSet;
  typedef std::vector<float> vfloat;
  typedef std::vector<std::string> vstring;
  typedef edm::View<reco::LeafCandidate> CandView;

  std::vector<edm::EDGetTokenT<edm::MergeableCounter> > eventCounterTokens_;

  FlatConsumers<bool> boolCSet_;
  FlatConsumers<int> intCSet_;
  FlatConsumers<double> doubleCSet_;
  FlatConsumers<float> floatCSet_;
  //FlatConsumers<std::string> stringCSet_;
  VectorConsumers<bool> vboolCSet_;
  VectorConsumers<int> vintCSet_;
  VectorConsumers<double> vdoubleCSet_;
  VectorConsumers<float> vfloatCSet_;
  VectorConsumers<std::string> vstringCSet_;

  CandConsumers candCSet_;

  TH1F* hNEvent_;

  TTree* tree_;
  int runNumber_, lumiNumber_, eventNumber_;

  struct FAILUREMODE
  {
    enum { KEEP, SKIP, ERROR };
  };
  int failureMode_;
};

GenericNtupleMaker::GenericNtupleMaker(const edm::ParameterSet& pset)
{
  std::string failureMode = pset.getUntrackedParameter<std::string>("failureMode", "keep");
  std::transform(failureMode.begin(), failureMode.end(), failureMode.begin(), ::tolower);
  if ( failureMode == "keep" ) failureMode_ = FAILUREMODE::KEEP;
  else if ( failureMode == "skip" ) failureMode_ = FAILUREMODE::SKIP;
  else if ( failureMode == "error" ) failureMode_ = FAILUREMODE::ERROR;
  else throw cms::Exception("ConfigError") << "select one from \"keep\", \"skip\", \"error\"\n";

  // Output histograms and tree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("event", "event");

  tree_->Branch("run"  , &runNumber_  , "run/I"  );
  tree_->Branch("lumi" , &lumiNumber_ , "lumi/I" );
  tree_->Branch("event", &eventNumber_, "event/I");

  boolCSet_.init(pset, "bool", consumesCollector(), tree_, "O");
  intCSet_.init(pset, "int", consumesCollector(), tree_, "I");
  doubleCSet_.init(pset, "double", consumesCollector(), tree_, "D");
  floatCSet_.init(pset, "float", consumesCollector(), tree_, "F");
  //stringCSet_.init(pset, "string", consumesCollector(), tree_, "F");
  vboolCSet_.init(pset, "bools", consumesCollector(), tree_);
  vintCSet_.init(pset, "ints", consumesCollector(), tree_);
  vdoubleCSet_.init(pset, "doubles", consumesCollector(), tree_);
  vfloatCSet_.init(pset, "floats", consumesCollector(), tree_);
  vstringCSet_.init(pset, "strings", consumesCollector(), tree_);

  candCSet_.init(pset, "cands", consumesCollector(), tree_);

  const auto eventCounters = pset.getParameter<vstring>("eventCounters");
  const size_t nEventCounter = eventCounters.size();
  hNEvent_ = fs->make<TH1F>("hNEvent", "NEvent", nEventCounter, 0, nEventCounter);
  for ( size_t i=0; i<nEventCounter; ++i )
  {
    hNEvent_->GetXaxis()->SetBinLabel(i+1, eventCounters[i].c_str());
    eventCounterTokens_.push_back(consumes<edm::MergeableCounter, edm::InLumi>(edm::InputTag(eventCounters[i])));
  }

}

void GenericNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  typedef edm::View<reco::LeafCandidate> Cands;

  int nFailure = 0;

  runNumber_   = event.run();
  lumiNumber_  = event.luminosityBlock();
  eventNumber_ = event.id().event();

  const bool doException = (failureMode_ == FAILUREMODE::ERROR);
  nFailure += boolCSet_.load(event, doException);
  nFailure += intCSet_.load(event, doException);
  nFailure += doubleCSet_.load(event, doException);
  nFailure += floatCSet_.load(event, doException);
  //nFailure += stringCSet_.load(event, doException);
  nFailure += vboolCSet_.load(event, doException);
  nFailure += vintCSet_.load(event, doException);
  nFailure += vdoubleCSet_.load(event, doException);
  nFailure += vfloatCSet_.load(event, doException);
  nFailure += vstringCSet_.load(event, doException);

  nFailure += candCSet_.load(event, doException);

  if ( nFailure == 0 or failureMode_ == FAILUREMODE::KEEP ) tree_->Fill();
  else if ( failureMode_ == FAILUREMODE::ERROR ) {
    edm::LogError("GenericNtupleMaker") << "Failed to get " << nFailure << " items";
    throw cms::Exception("DataError") << "Cannot get object from data";
  }
  //else if ( failureMode_ == FAILUREMODE::SKIP ); // don't fill and continue memory cleanup

  // Clear up after filling tree
  vboolCSet_.clear();
  vintCSet_.clear();
  vdoubleCSet_.clear();
  vfloatCSet_.clear();

  candCSet_.clear();

}

void GenericNtupleMaker::endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
  for ( size_t i=0, n=eventCounterTokens_.size(); i<n; ++i ) {
    edm::Handle<edm::MergeableCounter> eventCounterHandle;
    if ( lumi.getByToken(eventCounterTokens_[i], eventCounterHandle) ) {
      hNEvent_->Fill(i, double(eventCounterHandle->value));
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenericNtupleMaker);

