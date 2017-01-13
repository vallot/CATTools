import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
globalTag_rd   = '80X_dataRun2_2016SeptRepro_v4'
lumiJSON       = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON'
pileupMCmap    = '2016_25ns_Moriond17MC'

JetEnergyCorrection = ('Spring16_23Sep2016AllV2_DATA', 'Spring16_23Sep2016V2_MC')
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_DATA_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

