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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;
using namespace edm;

template<typename T>
class VectorConsumers
{
public:
  typedef edm::ParameterSet PSet;

  void init(const edm::ParameterSet& gpset, const string psetName, edm::ConsumesCollector && iC, TTree* tree)
  {
    if ( !gpset.existsAs<PSet>(psetName) ) return;
    const PSet pset = gpset.getParameter<PSet>(psetName);
    const auto names = pset.getParameterNamesForType<PSet>();
    for ( auto& name : names )
    {
      const auto ipset = pset.getParameter<PSet>(name);
      tokens_.push_back(iC.consumes<std::vector<T> >(ipset.getParameter<edm::InputTag>("src")));
      values_.push_back(new std::vector<T>);
      tree->Branch(name.c_str(), values_.back());
    }
  }

  int load(const edm::Event& event)
  {
    int nFailure = 0;
    for ( size_t i=0, n=tokens_.size(); i<n; ++i )
    {
      edm::Handle<std::vector<T> > handle;
      event.getByToken(tokens_[i], handle);
      if ( handle.isValid() )
      {
        values_[i]->insert(values_[i]->begin(), handle->begin(), handle->end());
      }
      else
      {
        ++nFailure;
      }
    }

    return nFailure;
  }

  void clear()
  {
    for ( auto& v : values_ ) v->clear();
  }

private:
  std::vector<edm::EDGetTokenT<std::vector<T> > > tokens_;
  std::vector<std::vector<T>*> values_;

};

template<typename T>
class FlatConsumers
{
public:
  typedef edm::ParameterSet PSet;

  void init(const edm::ParameterSet& gpset, const string psetName, edm::ConsumesCollector && iC,
            TTree* tree, const char* typeNameStr)
  {
    if ( !gpset.existsAs<PSet>(psetName) ) return;
    const PSet pset = gpset.getParameter<PSet>(psetName);
    const auto names = pset.getParameterNamesForType<PSet>();
    for ( auto& name : names )
    {
      const auto ipset = pset.getParameter<PSet>(name);
      tokens_.push_back(iC.consumes<T>(ipset.getParameter<edm::InputTag>("src")));
      values_.push_back(new T);
      tree->Branch(name.c_str(), values_.back(), (name+"/"+typeNameStr).c_str());
    }
  }

  int load(const edm::Event& event)
  {
    int nFailure = 0;
    for ( size_t i=0, n=tokens_.size(); i<n; ++i )
    {
      edm::Handle<T> handle;
      event.getByToken(tokens_[i], handle);
      if ( handle.isValid() ) *values_[i] = *handle;
      else
      {
        *values_[i] = -999;
        ++nFailure;
      }
    }

    return nFailure;
  }

private:
  std::vector<edm::EDGetTokenT<T> > tokens_;
  std::vector<T*> values_;

};

class GenericNtupleMaker : public edm::EDAnalyzer
{
public:
  GenericNtupleMaker(const edm::ParameterSet& pset);

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup) override;

private:
  typedef edm::ParameterSet PSet;
  typedef std::vector<int> vint;
  typedef std::vector<double> vdouble;
  typedef std::vector<std::string> strings;
  typedef edm::View<reco::LeafCandidate> CandView;
  typedef edm::ValueMap<double> Vmap;
  typedef edm::EDGetTokenT<CandView> CandToken;
  typedef edm::EDGetTokenT<Vmap> VmapToken;

  std::vector<CandToken> candTokens_;
  std::vector<std::vector<VmapToken> > vmapTokens_;
  std::vector<edm::EDGetTokenT<edm::MergeableCounter> > eventCounterTokens_;

  FlatConsumers<int> intCSet_;
  FlatConsumers<double> doubleCSet_;
  FlatConsumers<float> floatCSet_;
  VectorConsumers<int> vintCSet_;
  VectorConsumers<double> vdoubleCSet_;
  VectorConsumers<float> vfloatCSet_;

  typedef StringObjectFunction<reco::Candidate,true> CandFtn;
  typedef StringCutObjectSelector<reco::Candidate,true> CandSel;

  std::vector<int> indices_;
  std::vector<std::vector<CandFtn> > exprs_;
  std::vector<std::vector<CandSel> > selectors_;

  TH1F* hNEvent_;

  TTree* tree_;
  int runNumber_, lumiNumber_, eventNumber_;
  std::vector<std::vector<vdouble*> > candVars_;

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

  intCSet_.init(pset, "int", consumesCollector(), tree_, "I");
  doubleCSet_.init(pset, "double", consumesCollector(), tree_, "D");
  floatCSet_.init(pset, "float", consumesCollector(), tree_, "F");
  vintCSet_.init(pset, "ints", consumesCollector(), tree_);
  vdoubleCSet_.init(pset, "doubles", consumesCollector(), tree_);
  vfloatCSet_.init(pset, "floats", consumesCollector(), tree_);

  PSet candPSets = pset.getParameter<PSet>("cands");
  const strings candNames = candPSets.getParameterNamesForType<PSet>();
  for ( auto& candName : candNames )
  {
    PSet candPSet = candPSets.getParameter<PSet>(candName);

    edm::InputTag candToken = candPSet.getParameter<edm::InputTag>("src");
    candTokens_.push_back(consumes<CandView>(candToken));
    exprs_.push_back(std::vector<CandFtn>());
    selectors_.push_back(std::vector<CandSel>());
    vmapTokens_.push_back(std::vector<VmapToken>());
    candVars_.push_back(std::vector<vdouble*>());
    const string candTokenName = candToken.label();
    indices_.push_back(candPSet.getUntrackedParameter<int>("index", -1));
    const PSet exprSets = candPSet.getUntrackedParameter<PSet>("exprs", PSet());
    for ( auto& exprName : exprSets.getParameterNamesForType<string>() )
    {
      const string expr = exprSets.getParameter<string>(exprName);
      candVars_.back().push_back(new vdouble);
      exprs_.back().push_back(CandFtn(expr));

      tree_->Branch((candName+"_"+exprName).c_str(), candVars_.back().back());
    }
    const PSet selectionSets = candPSet.getUntrackedParameter<PSet>("seletions", PSet());
    for ( auto& selectionName : selectionSets.getParameterNamesForType<string>() )
    {
      const string selection = selectionSets.getParameter<string>(selectionName);
      candVars_.back().push_back(new vdouble);
      selectors_.back().push_back(CandSel(selection));

      tree_->Branch((candName+"_"+selectionName).c_str(), candVars_.back().back());
    }
    const strings vmapNames = candPSet.getUntrackedParameter<strings>("vmaps", strings());
    for ( auto& vmapName : vmapNames )
    {
      candVars_.back().push_back(new vdouble);

      edm::InputTag vmapToken(candTokenName, vmapName);
      vmapTokens_.back().push_back(consumes<Vmap>(vmapToken));

      tree_->Branch((candName+"_"+vmapName).c_str(), candVars_.back().back());
    }
  }

  const strings eventCounters = pset.getParameter<strings>("eventCounters");
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

  nFailure += intCSet_.load(event);
  nFailure += doubleCSet_.load(event);
  nFailure += floatCSet_.load(event);
  nFailure += vintCSet_.load(event);
  nFailure += vdoubleCSet_.load(event);
  nFailure += vfloatCSet_.load(event);

  const size_t nCand = candTokens_.size();
  for ( size_t iCand=0; iCand < nCand; ++iCand )
  {
    edm::Handle<CandView> srcHandle;
    event.getByToken(candTokens_[iCand], srcHandle);
    if ( !srcHandle.isValid() )
    {
      ++nFailure;
      continue;
    }

    const int index = indices_[iCand];
    const std::vector<CandFtn>& exprs = exprs_[iCand];
    const std::vector<CandSel>& selectors = selectors_[iCand];
    std::vector<VmapToken>& vmapTokens = vmapTokens_[iCand];
    const size_t nExpr = exprs.size();
    const size_t nSels = selectors.size();
    const size_t nVmap = vmapTokens.size();
    std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVmap);
    for ( size_t iVar=0; iVar<nVmap; ++iVar )
    {
      event.getByToken(vmapTokens[iVar], vmapHandles[iVar]);
    }

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      if ( index >= 0 and int(i) != index ) continue;
      edm::Ref<CandView> candRef(srcHandle, i);

      for ( size_t j=0; j<nExpr; ++j )
      {
        const double val = exprs[j](*candRef);
        candVars_[iCand][j]->push_back(val);
      }
      for ( size_t j=0; j<nSels; ++j )
      {
        const double val = selectors[j](*candRef);
        candVars_[iCand][j+nExpr]->push_back(val);
      }
      for ( size_t j=0; j<nVmap; ++j )
      {
        double val = 0;
        if ( vmapHandles[j].isValid() ) val = (*vmapHandles[j])[candRef];
        candVars_[iCand][j+nExpr+nSels]->push_back(val);
      }
    }
  }

  if ( nFailure == 0 or failureMode_ == FAILUREMODE::KEEP ) tree_->Fill();
  else if ( failureMode_ == FAILUREMODE::ERROR )
  {
    edm::LogError("GenericNtupleMaker") << "Failed to get " << nFailure << " items";
    throw cms::Exception("DataError") << "Cannot get object from data";
  }
  //else if ( failureMode_ == FAILUREMODE::SKIP ); // don't fill and continue memory cleanup

  // Clear up after filling tree
  vintCSet_.clear();
  vdoubleCSet_.clear();

  for ( size_t iCand=0; iCand<nCand; ++iCand )
  {
    const size_t nVar = candVars_[iCand].size();
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      candVars_[iCand][iVar]->clear();
    }
  }
}

void GenericNtupleMaker::endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
  for ( size_t i=0, n=eventCounterTokens_.size(); i<n; ++i )
  {
    edm::Handle<edm::MergeableCounter> eventCounterHandle;
    if ( lumi.getByToken(eventCounterTokens_[i], eventCounterHandle) )
    {
      hNEvent_->Fill(i, double(eventCounterHandle->value));
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenericNtupleMaker);

