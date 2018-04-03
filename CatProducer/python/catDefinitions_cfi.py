import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
globalTag_rd   = '80X_dataRun2_2016SeptRepro_v7'
lumiJSON       = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON'
pileupMCmap    = 'mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU'

JetEnergyCorrection = ('Summer16_23Sep2016AllV4_DATA', 'Summer16_23Sep2016V4_MC')
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection[1]
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

