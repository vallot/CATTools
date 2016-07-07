#ifndef __CATTools_TopDileptonEventCommon__
#define __CATTools_TopDileptonEventCommon__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TH1D.h"

#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"
#include "CATTools/CatAnalyzer/interface/TopEventCommon.h"
class TopDileptonEventCommon : public TopEventCommon {
public:
  explicit TopDileptonEventCommon(const edm::ParameterSet&);
  virtual ~TopDileptonEventCommon(){ showSummary(); }
  virtual void analyzeCustom(const edm::Event&, const edm::EventSetup&, int sys );
  virtual void showSummary();

// Use protect keyword for branches.
protected : 

private:

};


#endif
