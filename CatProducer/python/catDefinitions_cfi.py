import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
#globalTag_mc   = '103X_upgrade2018_realistic_v8'
#globalTag_rd   = '103X_dataRun2_Prompt_v3'
globalTag_mc   = '102X_upgrade2018_realistic_v12'
globalTag_rd   = '102X_dataRun2_Sep2018Rereco_v1'
#globalTag_rd   = '102X_dataRun2_Prompt_v11'
lumiJSON       = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON'
pileupMCmap    = '2018_25ns_MC'

JetEnergyCorrection = ('Fall17_17Nov2017_V32_94X_DATA', 'Fall17_17Nov2017_V32_94X_MC')
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection[1]
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

