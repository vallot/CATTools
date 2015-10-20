from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'PAT2CAT_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.publishDataName = 'cat'

config.section_("Site")
# Where the output files will be transmitted to
#config.Site.storageSite = 'T2_KR_KNU'
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/CAT/
#config.Site.storageSite = 'T3_KR_KISTI'
config.Site.storageSite = 'T3_KR_UOS'
config.Data.outLFNDirBase = '/store/group/CAT/' 
#config.Site.storageSite = 'T3_US_FNALLPC'
