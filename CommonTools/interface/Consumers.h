#ifndef CATTools_CommonTools_Consumers_H
#define CATTools_CommonTools_Consumers_H

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TTree.h"

#include <memory>
#include <vector>
#include <string>

namespace cat {

template<typename T>
class VectorConsumers
{
public:
  typedef edm::ParameterSet PSet;

  void init(const edm::ParameterSet& gpset, const std::string psetName, edm::ConsumesCollector && iC, TTree* tree)
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
    const auto labels = pset.getParameterNamesForType<edm::InputTag>();
    for ( auto& name : labels )
    {
      tokens_.push_back(iC.consumes<std::vector<T> >(pset.getParameter<edm::InputTag>(name)));
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

  void init(const edm::ParameterSet& gpset, const std::string psetName, edm::ConsumesCollector && iC,
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
    const auto labels = pset.getParameterNamesForType<edm::InputTag>();
    for ( auto& name : labels )
    {
      tokens_.push_back(iC.consumes<T>(pset.getParameter<edm::InputTag>(name)));
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

}

#endif
