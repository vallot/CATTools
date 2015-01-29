from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'runCat_cfg.py'
config.JobType.pyCfgParams = ['runOnMC=False','globalTag=FT53_V21A_AN6::All']
config.JobType.inputFiles  = ['Winter14_V5_DATA.db']

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 20
config.Data.lumiMask     = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'
config.Data.publication  = config.General.transferOutputs
config.Data.publishDataName = 'cat'

config.section_("Site")
# Where the output files will be transmitted to
#config.Site.storageSite = 'T2_KR_KNU'
config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T3_KR_UOS'
#config.Site.storageSite = 'T3_US_FNALLPC'

