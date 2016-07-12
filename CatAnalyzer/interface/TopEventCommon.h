#ifndef __CATTools_TopEventCommon__
#define __CATTools_TopEventCommon__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TH1D.h"

#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"

class TopEventCommon : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit TopEventCommon(const edm::ParameterSet&);
  virtual ~TopEventCommon();
  void paramInit(const edm::ParameterSet&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void analyzeCustom(const edm::Event&, const edm::EventSetup&, int sys );
  void genInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup );
  void setBranch(TTree* tree, int sys);
  virtual void setBranchCustom(TTree* tree, int sys) {}
  virtual void resetBranchCustom() {}
  void resetBr() {
    TopEventInfo& evInfo_ = TopEventInfo::getInstance();
    evInfo_.resetBr();
    resetBranchCustom();
  }
  int runEventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys) {
   return eventSelect_->eventSelection(iEvent, iSetup, tree, sys);
  }
  virtual void setEventSelection( const edm::ParameterSet& iConfig  ) { 
    eventSelect_ = new TTEventSelector( iConfig, consumesCollector());
  }
  virtual void showSummary();
  virtual void beginJob();
  virtual void endJob();
// Use protect keyword for branches.
protected : 
  TTEventSelector* eventSelect_;
  // I/O variables.
  std::vector<TTree*> ttree_;
  TH1D * h_nevents;
  const edm::ParameterSet& iConfig_;

  // Token for gen Information
  edm::EDGetTokenT<float> genWeightToken_;
  edm::EDGetTokenT<std::vector<float>> pdfweightToken_, scaleupweightsToken_, scaledownweightsToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightToken_up_, puweightToken_dn_, topPtWeight_;

  edm::EDGetTokenT<int>          partonTop_channel_;
  edm::EDGetTokenT<std::vector<int> > partonTop_modes_;
  edm::EDGetTokenT<reco::GenParticleCollection> partonTop_genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > pseudoTop_leptons_, pseudoTop_neutrinos_, pseudoTop_jets_;

private:
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&);
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};


};


#endif
