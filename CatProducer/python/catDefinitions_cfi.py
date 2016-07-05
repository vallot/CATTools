import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '80X_mcRun2_asymptotic_v14'
globalTag_rd   = '80X_dataRun2_Prompt_v9'
lumiJSON       = 'Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON_NoL1T'
pileupMCmap    = '2016_25ns_SpringMC'

JetEnergyCorrection = 'Spring16_25nsV6'
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_DATA_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_DATA_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection

