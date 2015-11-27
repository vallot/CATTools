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
  for ( auto c : candVars_ )
  {
    for ( auto v : c )
    {
      delete v;
    }
  }
}

void CandConsumers::init(const edm::ParameterSet& gpset, const std::string psetName, edm::ConsumesCollector&& iC, TTree* tree)
{
  if ( !gpset.existsAs<PSet>(psetName) ) return;
  const auto pset = gpset.getParameter<PSet>(psetName);
  const auto candNames = pset.getParameterNamesForType<PSet>();
  for ( auto& candName : candNames )
  {
    PSet candPSet = pset.getParameter<PSet>(candName);

    edm::InputTag candToken = candPSet.getParameter<edm::InputTag>("src");
    candTokens_.push_back(iC.consumes<CandView>(candToken));
    exprs_.push_back(std::vector<CandFtn>());
    boolexprs_.push_back(std::vector<CandSel>());
    vmapTokens_.push_back(std::vector<VmapToken>());
    candVars_.push_back(std::vector<vfloat*>());
    const string candTokenName = candToken.label();
    indices_.push_back(candPSet.getUntrackedParameter<int>("index", -1));
    const PSet exprSets = candPSet.getUntrackedParameter<PSet>("exprs", PSet());
    for ( auto& exprName : exprSets.getParameterNamesForType<string>() )
    {
      const string expr = exprSets.getParameter<string>(exprName);
      candVars_.back().push_back(new vfloat);
      exprs_.back().push_back(CandFtn(expr));

      if ( tree ) tree->Branch((candName+"_"+exprName).c_str(), candVars_.back().back());
    }
    const PSet boolexprSets = candPSet.getUntrackedParameter<PSet>("boolexprs", PSet());
    for ( auto& boolexprName : boolexprSets.getParameterNamesForType<string>() )
    {
      const string boolexpr = boolexprSets.getParameter<string>(boolexprName);
      candVars_.back().push_back(new vfloat);
      boolexprs_.back().push_back(CandSel(boolexpr));

      if ( tree ) tree->Branch((candName+"_"+boolexprName).c_str(), candVars_.back().back());
    }
    const auto vmapNames = candPSet.getUntrackedParameter<vstring>("vmaps", vstring());
    for ( auto& vmapName : vmapNames )
    {
      candVars_.back().push_back(new vfloat);

      edm::InputTag vmapToken(candTokenName, vmapName);
      vmapTokens_.back().push_back(iC.consumes<Vmap>(vmapToken));

      if ( tree ) tree->Branch((candName+"_"+vmapName).c_str(), candVars_.back().back());
    }
  }
}

int CandConsumers::load(const edm::Event& event, const bool doException)
{
  int nFailure = 0;
  const size_t nCand = candTokens_.size();
  for ( size_t iCand=0; iCand < nCand; ++iCand )
  {
    edm::Handle<CandView> srcHandle;
    event.getByToken(candTokens_[iCand], srcHandle);
    if ( !srcHandle.isValid() )
    {
      if ( doException ) throw cms::Exception("DataError") << "Cannot load " << srcHandle;
      ++nFailure;
      continue;
    }

    const int index = indices_[iCand];
    const std::vector<CandFtn>& exprs = exprs_[iCand];
    const std::vector<CandSel>& boolexprs = boolexprs_[iCand];
    std::vector<VmapToken>& vmapTokens = vmapTokens_[iCand];
    const size_t nExpr = exprs.size();
    const size_t nBExprs = boolexprs.size();
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
        double val = exprs[j](*candRef);
        if ( !std::isfinite(val) ) val = -999;
        candVars_[iCand][j]->push_back(val);
      }
      for ( size_t j=0; j<nBExprs; ++j )
      {
        const double val = boolexprs[j](*candRef);
        candVars_[iCand][j+nExpr]->push_back(val);
      }
      for ( size_t j=0; j<nVmap; ++j )
      {
        double val = 0;
        if ( vmapHandles[j].isValid() ) val = (*vmapHandles[j])[candRef];
        else if ( doException ) throw cms::Exception("DataError") << "Cannot load " << vmapHandles[j];
        candVars_[iCand][j+nExpr+nBExprs]->push_back(val);
      }
    }
  }

  return nFailure;
}

void CandConsumers::clear()
{
  const size_t nCand = candVars_.size();
  for ( size_t iCand=0; iCand<nCand; ++iCand )
  {
    const size_t nVar = candVars_[iCand].size();
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      candVars_[iCand][iVar]->clear();
    }
  }
}



}
