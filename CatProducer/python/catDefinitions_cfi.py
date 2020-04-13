import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '102X_mcRun2_asymptotic_v7'
globalTag_rd   = '102X_dataRun2_v12'
lumiJSON       = 'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON'
pileupMCmap    = '2016_25ns_Moriond17MC'

JetEnergyCorrection = ('Summer16_07Aug2017All_V11_DATA', 'Summer16_07Aug2017_V11_MC')
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection[1]
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

