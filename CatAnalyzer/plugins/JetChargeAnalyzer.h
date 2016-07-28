#ifndef __JetChargeAnalyzer__
#define __JetChargeAnalyzer__

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
//#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include<TTree.h>
#include<TNtuple.h>

struct Data {
  // Variables.
  float lep_pt[2];
  int lep_pdgId[2];
  
  float jet_pt[2];
  int jet_pdgId[2];
  int jet_charge[2];
  Data(){
    reset();
  }
  Data(TTree* tree){
    reset();
    bookingTree( tree ); 
  }
  void reset() {
    for( int i= 0; i< 2 ; ++i) {
      lep_pt[i]= 0.0f;
      lep_pdgId[i]= 0;
      jet_pt[i] = 0.0f;
      jet_pdgId[i]= 0;
      jet_charge[i]= 0;
    }
  }

  void bookingTree( TTree* tree) {
    tree->Branch("data",&(this->lep_pt[0]),"lep_pt[2]/F:lep_pdgId[2]/I:jet_pt[2]/F:jet_pdgId[2]/I:jet_charge[2]/I");
  }
};

class JetChargeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit JetChargeAnalyzer(const edm::ParameterSet&) ;
  ~JetChargeAnalyzer(){};
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons) const;
  float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs) const;

  cat::JetCollection selectJets(const cat::JetCollection& jets, const std::vector<const cat::Lepton*>& recolep);
  cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
// Use protect keyword for branches.
protected : 
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::MuonCollection>      muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      elecToken_;
  TH1D* h_nevents;
  TTree* ttree_;
  Data* data;
  int b_run;
  int b_event;

private :
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override {};
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};
  //ScaleFactorEvaluator muonSF_, elecSF_;


};


#endif
