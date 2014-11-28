from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'cat_test1'
config.General.transferLogs = False
config.General.transferOutputs = False

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'runCat.py'
config.JobType.pyCfgParams = ['runOnMC=True']

config.section_("Data")
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset = '/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 5000
config.Data.publication = False

# This string is used to construct the output dataset name
config.Data.publishDataName = 'cat'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_KR_KNU'
#config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T3_KR_UOS'
#config.Site.storageSite = 'T3_US_FNALLPC'


