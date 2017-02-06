import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
globalTag_rd   = '80X_dataRun2_2016SeptRepro_v7'
lumiJSON       = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON'
pileupMCmap    = '2016_25ns_Moriond17MC'

JetEnergyCorrection = ('Summer16_23Sep2016AllV4_DATA', 'Summer16_23Sep2016V4_MC')
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection[1]
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

