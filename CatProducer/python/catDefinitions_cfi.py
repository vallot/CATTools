import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '102X_upgrade2018_realistic_v20'
globalTag_rd   = '102X_dataRun2_v12'
#globalTag_rd   = '102X_dataRun2_Prompt_v15'
lumiJSON       = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON'
pileupMCmap    = '2018_25ns_MC'

JetEnergyCorrection = ('Autumn18_RunABCD_V19_DATA', 'Autumn18_V19_MC')
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection[1]
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

