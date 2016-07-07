#ifndef __CATTools_CatAnalyzer_TTDileptonEventSelector__
#define __CATTools_CatAnalyzer_TTDileptonEventSelector__

#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"

class TTDileptonEventSelector : public TTEventSelector 
{
public :
  explicit TTDileptonEventSelector(const edm::ParameterSet&, edm::ConsumesCollector&& iC );
  ~TTDileptonEventSelector() {}
  virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys) ;
  void setBranch(TTree* tree, int sys);
  void resetBranch(); 

  /*
  virtual float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, TopEventCommonGlobal::sys_e sys) const;
  virtual float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, TopEventCommonGlobal::sys_e sys) const;
  virtual cat::JetCollection selectJets(const cat::JetCollection& jets, const TopEventCommonGlobal::LeptonPtrs& recolep, TopEventCommonGlobal::sys_e sys);
  virtual cat::JetCollection selectBJets(const cat::JetCollection& jets) const;
  */

protected:
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_;

};
#endif
