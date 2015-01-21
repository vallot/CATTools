#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
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

class GenericNtupleMaker : public edm::EDAnalyzer
{
public:
  GenericNtupleMaker(const edm::ParameterSet& pset);

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup) override;

private:
  typedef edm::ParameterSet PSet;
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;
  typedef edm::View<reco::LeafCandidate> CandView;
  typedef edm::ValueMap<double> Vmap;
  typedef edm::EDGetTokenT<CandView> CandToken;
  typedef edm::EDGetTokenT<Vmap> VmapToken;

  std::vector<CandToken> candTokens_;
  std::vector<std::vector<VmapToken> > vmapTokens_;
  std::vector<edm::EDGetTokenT<int> > intTokens_;
  std::vector<edm::EDGetTokenT<double> > doubleTokens_;
  std::vector<edm::EDGetTokenT<doubles> > doublesTokens_;
  std::vector<edm::EDGetTokenT<edm::MergeableCounter> > eventCounterTokens_;

  typedef StringObjectFunction<reco::Candidate,true> CandFtn;
  typedef StringCutObjectSelector<reco::Candidate,true> CandSel;

  std::vector<int> indices_;
  std::vector<std::vector<CandFtn> > exprs_;
  std::vector<std::vector<CandSel> > selectors_;

  TH1F* hNEvent_;

  TTree* tree_;
  int runNumber_, lumiNumber_, eventNumber_;
  std::vector<int*> ints_;
  std::vector<double*> double_;
  std::vector<doubles*> vdouble_;
  std::vector<std::vector<doubles*> > candVars_;

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

  PSet intPSets = pset.getParameter<PSet>("int");
  const strings intNames = intPSets.getParameterNamesForType<PSet>();
  for ( auto& intName : intNames )
  {
    PSet intPSet = intPSets.getParameter<PSet>(intName);
    intTokens_.push_back(consumes<int>(intPSet.getParameter<edm::InputTag>("src")));

    ints_.push_back(new int);
    tree_->Branch(intName.c_str(), ints_.back(), (intName+"/I").c_str());
  }

  PSet doublePSets = pset.getParameter<PSet>("double");
  const strings doubleNames = doublePSets.getParameterNamesForType<PSet>();
  for ( auto& doubleName : doubleNames )
  {
    PSet doublePSet = doublePSets.getParameter<PSet>(doubleName);
    doubleTokens_.push_back(consumes<double>(doublePSet.getParameter<edm::InputTag>("src")));

    double_.push_back(new double);
    tree_->Branch(doubleName.c_str(), double_.back(), (doubleName+"/D").c_str());
  }

  PSet doublesPSets = pset.getParameter<PSet>("doubles");
  const strings doublesNames = doublesPSets.getParameterNamesForType<PSet>();
  for ( auto& doublesName : doublesNames )
  {
    PSet doublesPSet = doublesPSets.getParameter<PSet>(doublesName);
    doublesTokens_.push_back(consumes<doubles>(doublesPSet.getParameter<edm::InputTag>("src")));

    vdouble_.push_back(new doubles);
    tree_->Branch(doublesName.c_str(), vdouble_.back());
  }

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
    candVars_.push_back(std::vector<doubles*>());
    const string candTokenName = candToken.label();
    indices_.push_back(candPSet.getUntrackedParameter<int>("index", -1));
    const PSet exprSets = candPSet.getUntrackedParameter<PSet>("exprs", PSet());
    for ( auto& exprName : exprSets.getParameterNamesForType<string>() )
    {
      const string expr = exprSets.getParameter<string>(exprName);
      candVars_.back().push_back(new doubles);
      exprs_.back().push_back(CandFtn(expr));

      tree_->Branch((candName+"_"+exprName).c_str(), candVars_.back().back());
    }
    const PSet selectionSets = candPSet.getUntrackedParameter<PSet>("seletions", PSet());
    for ( auto& selectionName : selectionSets.getParameterNamesForType<string>() )
    {
      const string selection = selectionSets.getParameter<string>(selectionName);
      candVars_.back().push_back(new doubles);
      selectors_.back().push_back(CandSel(selection));

      tree_->Branch((candName+"_"+selectionName).c_str(), candVars_.back().back());
    }
    const strings vmapNames = candPSet.getUntrackedParameter<strings>("vmaps", strings());
    for ( auto& vmapName : vmapNames )
    {
      candVars_.back().push_back(new doubles);

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

  for ( size_t i=0, n=intTokens_.size(); i<n; ++i )
  {
    edm::Handle<int> intHandle;
    event.getByToken(intTokens_[i], intHandle);
    if ( intHandle.isValid() ) *ints_[i] = *intHandle;
    else
    {
      *ints_[i] = 0;
      ++nFailure;
    }
  }

  for ( size_t i=0, n=doubleTokens_.size(); i<n; ++i )
  {
    edm::Handle<double> doubleHandle;
    event.getByToken(doubleTokens_[i], doubleHandle);

    if ( doubleHandle.isValid() ) *double_[i] = *doubleHandle;
    else
    {
      *double_[i] = 0;
      ++nFailure;
    }
  }

  for ( size_t i=0, n=doublesTokens_.size(); i<n; ++i )
  {
    edm::Handle<doubles> doublesHandle;
    event.getByToken(doublesTokens_[i], doublesHandle);

    if ( doublesHandle.isValid() )
    {
      vdouble_[i]->insert(vdouble_[i]->begin(), doublesHandle->begin(), doublesHandle->end());
    }
    else
    {
      ++nFailure;
    }
  }

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

  for ( size_t i=0, n=doublesTokens_.size(); i<n; ++i )
  {
    edm::Handle<doubles> doublesHandle;
    event.getByToken(doublesTokens_[i], doublesHandle);

    vdouble_[i]->clear();
  }

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


DEFINE_FWK_MODULE(GenericNtupleMaker);

