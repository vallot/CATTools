#ifndef __CATTools_CatAnalyzer_TTDiLeptonEventSelector__
#define __CATTools_CatAnalyzer_TTDiLeptonEventSelector__

#include "CATTools/CatAnalyzer/interface/TTEventSelector.h"

class TTDiLeptonEventSelector : public TTEventSelector 
{
public :
  explicit TTDiLeptonEventSelector(const edm::ParameterSet&, edm::ConsumesCollector&& iC );
  ~TTDiLeptonEventSelector() {}
  virtual int eventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup, TTree* tree, int sys) ;
  void setBranch(TTree* tree, int sys);
  void resetBranch(); 

  virtual float selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, TopEventCommonGlobal::sys_e sys) const;
  virtual float selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, TopEventCommonGlobal::sys_e sys) const;

protected:
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_;

};
#endif
