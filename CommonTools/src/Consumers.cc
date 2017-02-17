#include "CATTools/CommonTools/interface/Consumers.h"

using namespace std;

namespace cat {
// Build dummy instances here
FlatConsumers<bool> boolCSet;
FlatConsumers<int> intCSet;
FlatConsumers<double> doubleCSet;
FlatConsumers<float> floatCSet;
FlatConsumers<std::string> stringCSet;

VectorConsumers<bool> vboolCSet;
VectorConsumers<int> vintCSet;
VectorConsumers<double> vdoubleCSet;
VectorConsumers<float> vfloatCSet;
VectorConsumers<std::string> vstringCSet;

// CandConsumer implementation
CandConsumers::~CandConsumers()
{
  for ( auto c : candVarsF_ ) { for ( auto v : c ) delete[] v; }
  for ( auto c : candVarsI_ ) { for ( auto v : c ) delete[] v; }
  for ( auto c : candVarsB_ ) { for ( auto v : c ) delete[] v; }
  for ( auto c : candSize_ ) delete c;
}

void CandConsumers::init(const edm::ParameterSet& gpset, const std::string psetName, edm::ConsumesCollector&& iC, TTree* tree)
{
  if ( !tree ) return; // no need to do this step if tree is invalid
  if ( !gpset.existsAs<PSet>(psetName) ) return;
  const auto pset = gpset.getParameter<PSet>(psetName);
  const auto candNames = pset.getParameterNamesForType<PSet>();
  for ( auto& candName : candNames ) {
    PSet candPSet = pset.getParameter<PSet>(candName);

    edm::InputTag candToken = candPSet.getParameter<edm::InputTag>("src");
    candTokens_.push_back(iC.consumes<CandView>(candToken));
    exprsF_.push_back(std::vector<CandFtn>());
    exprsI_.push_back(std::vector<CandFtn>());
    exprsB_.push_back(std::vector<CandFtn>());
    boolexprs_.push_back(std::vector<CandSel>());
    vmapTokens_.push_back(std::vector<VmapToken>());

    candSize_.push_back(new unsigned short);
    candVarsF_.push_back(std::vector<float*>());
    candVarsI_.push_back(std::vector<int*>());
    candVarsB_.push_back(std::vector<bool*>());

    const string sizeName = candName+"_n";
    tree->Branch(sizeName.c_str(), candSize_.back(), (sizeName+"/s").c_str());

    const string candTokenName = candToken.label();
    indices_.push_back(candPSet.getUntrackedParameter<int>("index", -1));
    const PSet exprFSets = candPSet.getUntrackedParameter<PSet>("exprsF", PSet());
    for ( auto& exprName : exprFSets.getParameterNamesForType<string>() ) {
      const string expr = exprFSets.getParameter<string>(exprName);
      candVarsF_.back().push_back(new float[maxSize_]);
      exprsF_.back().push_back(CandFtn(expr));

      const string brName = candName+"_"+exprName;
      tree->Branch(brName.c_str(), candVarsF_.back().back(), (brName+"["+sizeName+"]/F").c_str());
    }
    const PSet exprISets = candPSet.getUntrackedParameter<PSet>("exprsI", PSet());
    for ( auto& exprName : exprISets.getParameterNamesForType<string>() ) {
      const string expr = exprISets.getParameter<string>(exprName);
      candVarsI_.back().push_back(new int[maxSize_]);
      exprsI_.back().push_back(CandFtn(expr));

      const string brName = candName+"_"+exprName;
      tree->Branch(brName.c_str(), candVarsI_.back().back(), (brName+"["+sizeName+"]/I").c_str());
    }
    const PSet exprBSets = candPSet.getUntrackedParameter<PSet>("exprsB", PSet());
    for ( auto& exprName : exprBSets.getParameterNamesForType<string>() ) {
      const string expr = exprBSets.getParameter<string>(exprName);
      candVarsB_.back().push_back(new bool[maxSize_]);
      exprsB_.back().push_back(CandFtn(expr));

      const string brName = candName+"_"+exprName;
      tree->Branch(brName.c_str(), candVarsB_.back().back(), (brName+"["+sizeName+"]/O").c_str());
    }
    const PSet boolexprSets = candPSet.getUntrackedParameter<PSet>("boolexprs", PSet());
    for ( auto& exprName : boolexprSets.getParameterNamesForType<string>() ) {
      const string expr = boolexprSets.getParameter<string>(exprName);
      candVarsB_.back().push_back(new bool[maxSize_]);
      boolexprs_.back().push_back(CandSel(expr));

      const string brName = candName+"_"+exprName;
      tree->Branch(candName.c_str(), candVarsB_.back().back(), (brName+"["+sizeName+"]/O").c_str());
    }
    const auto vmapNames = candPSet.getUntrackedParameter<vstring>("vmaps", vstring());
    for ( auto& vmapName : vmapNames ) {
      candVarsF_.back().push_back(new float[maxSize_]);

      edm::InputTag vmapToken(candTokenName, vmapName);
      vmapTokens_.back().push_back(iC.consumes<Vmap>(vmapToken));

      const string brName = candName+"_"+vmapName;
      tree->Branch(brName.c_str(), candVarsF_.back().back(), (brName+"["+sizeName+"]/F").c_str());
    }
  }
}

int CandConsumers::load(const edm::Event& event, const bool doException)
{
  int nFailure = 0;
  const size_t nCand = candTokens_.size();
  for ( size_t iCand=0; iCand < nCand; ++iCand ) {
    edm::Handle<CandView> srcHandle;
    event.getByToken(candTokens_[iCand], srcHandle);
    if ( !srcHandle.isValid() ) {
      if ( doException ) throw cms::Exception("DataError") << "Cannot load " << srcHandle;
      ++nFailure;
      continue;
    }

    const int index = indices_[iCand];
    const std::vector<CandFtn>& exprsF = exprsF_[iCand];
    const std::vector<CandFtn>& exprsI = exprsI_[iCand];
    const std::vector<CandFtn>& exprsB = exprsB_[iCand];
    const std::vector<CandSel>& boolexprs = boolexprs_[iCand];
    std::vector<VmapToken>& vmapTokens = vmapTokens_[iCand];
    const size_t nExprF = exprsF.size();
    const size_t nExprI = exprsI.size();
    const size_t nExprB = exprsB.size();
    const size_t nBExprs = boolexprs.size();
    const size_t nVmap = vmapTokens.size();
    std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVmap);
    for ( size_t iVar=0; iVar<nVmap; ++iVar ) {
      event.getByToken(vmapTokens[iVar], vmapHandles[iVar]);
    }

    int candSize = 0;
    for ( size_t i=0, n=srcHandle->size(); i<n and candSize < maxSize_; ++i ) {
      if ( index >= 0 and int(i) != index ) continue;
      edm::Ref<CandView> candRef(srcHandle, i);

      for ( size_t j=0; j<nExprF; ++j ) {
        double val = exprsF[j](*candRef);
        if ( !std::isfinite(val) ) val = -999;
        candVarsF_[iCand][j][candSize] = val;
      }
      for ( size_t j=0; j<nExprI; ++j ) {
        const int val = exprsI[j](*candRef);
        candVarsI_[iCand][j][candSize] = val;
      }
      for ( size_t j=0; j<nExprB; ++j ) {
        const bool val = exprsB[j](*candRef);
        candVarsB_[iCand][j][candSize] = val;
      }
      for ( size_t j=0; j<nBExprs; ++j ) {
        const bool val = boolexprs[j](*candRef);
        candVarsB_[iCand][j+nExprB][candSize] = val;
      }
      for ( size_t j=0; j<nVmap; ++j ) {
        double val = 0;
        if ( vmapHandles[j].isValid() ) val = (*vmapHandles[j])[candRef];
        else if ( doException ) throw cms::Exception("DataError") << "Cannot load " << vmapHandles[j];
        candVarsF_[iCand][j+nExprF][candSize] = val;
      }

      ++candSize;
    }
    *candSize_[iCand] = candSize;
  }

  return nFailure;
}

void CandConsumers::clear()
{
  for ( auto& x : candSize_ ) *x = 0;
}

}
