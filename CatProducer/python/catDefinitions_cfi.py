import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '76X_mcRun2_asymptotic_v12'
globalTag_rd   = '76X_dataRun2_v15'
lumiJSON       = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2'
lumiJSONSilver = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2'
pileupMCmap    = '2015_25ns_FallMC'

JetEnergyCorrection = 'Summer15_25nsV6'
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_DATA_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection

