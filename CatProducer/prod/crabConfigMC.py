from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'PAT2CAT_cfg.py'
config.JobType.pyCfgParams = ['runOnMC=True',"globalTag='PHYS14_25_V2::All'"]
config.JobType.inputFiles  = ['CATTools/CatProducer/data/PHYS14_V4_MC.db']

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 1
config.Data.totalUnits   = 5000
config.Data.publication  = True
config.Data.publishDataName = 'cat'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_KR_KNU'
#config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T3_KR_UOS'
#config.Site.storageSite = 'T3_US_FNALLPC'
