import FWCore.ParameterSet.Config as cms

bunchCrossing  = 25
globalTag_mc   = '94X_mc2017_realistic_v14'
globalTag_rd   = '94X_dataRun2_v6'
lumiJSON       = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON'
pileupMCmap    = '2017_25ns_WinterMC'

JetEnergyCorrection = ('Summer16_23Sep2016AllV4_DATA', 'Summer16_23Sep2016V4_MC')
JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_UncertaintySources_AK4PFchs.txt'%JetEnergyCorrection[1]
#JECUncertaintyFile  = 'CATTools/CatProducer/data/JEC/%s_Uncertainty_AK4PFchs.txt'%JetEnergyCorrection[1]

