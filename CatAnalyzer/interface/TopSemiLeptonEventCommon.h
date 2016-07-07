#ifndef __CATTools_TopSemiLeptonEventCommon__
#define __CATTools_TopSemiLeptonEventCommon__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TH1D.h"

#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"
#include "CATTools/CatAnalyzer/interface/TopEventCommon.h"

class TopSemiLeptonEventCommon : public TopEventCommon {
public:
  explicit TopSemiLeptonEventCommon(const edm::ParameterSet&);
  virtual ~TopSemiLeptonEventCommon() { showSummary(); }
  virtual void analyzeCustom(const edm::Event&, const edm::EventSetup&, int sys );
  virtual void showSummary();

// Use protect keyword for branches.
protected : 

private:

};


#endif
