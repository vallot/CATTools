#ifndef __CATTools_CatAnalyzer_TTSemiLeptonEventSelector__
#define __CATTools_CatAnalyzer_TTSemiLeptonEventSelector__

#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"

class TTSemiLeptonEventSelector : public TTEventSelector 
{
public :
  explicit TTSemiLeptonEventSelector(const edm::ParameterSet&, edm::ConsumesCollector&& iC );
  //explicit TTSemiLeptonEventSelector(const edm::ParameterSet&, edm::ConsumesCollector iC );
  ~TTSemiLeptonEventSelector() {}
  virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys) override;
  void setBranch(TTree* tree, int sys) override;
  void resetBranch() override; 

  virtual float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, cat::MuonCollection& vetomuons, TopEventCommonGlobal::sys_e sys) const ;
  virtual float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, cat::ElectronCollection& vetoelecs, TopEventCommonGlobal::sys_e sys) const ;
  /*
  virtual cat::JetCollection selectJets(const cat::JetCollection& jets, const TopEventCommonGlobal::LeptonPtrs& recolep, TopEventCommonGlobal::sys_e sys);
  virtual cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  */
  bool isVetoMuon(cat::Muon) const;
  bool isVetoElec(cat::Electron) const;
  bool isSelectMuon(cat::Muon) const;
  bool isSelectElec(cat::Electron) const;


protected:
  //edm::EDGetTokenT<std::string> vetoMuonIDCut_, vetoElectornIDCut_;
  std::string vetoElectronIDCut_ ;

  float vetoMuonPtCut_, vetoMuonEtaCut_, vetoMuonIsoCut_;
  float vetoElectronPtCut_, vetoElectronEtaCut_, vetoElectronIsoCut_;

  edm::EDGetTokenT<int> trigTokenMUJET_;
  edm::EDGetTokenT<int> trigTokenELJET_;

};
#endif
