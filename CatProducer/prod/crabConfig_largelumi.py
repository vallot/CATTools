from WMCore.Configuration import Configuration
from CRABAPI.RawCommand import crabCommand
from dbs.apis.dbsClient import DbsApi
dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')

dataset = '/W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v11-v1/MINIAODSIM'
fileDictList=dbs.listFiles(dataset=dataset)
lfnList = [ dic['logical_file_name'] for dic in fileDictList ]

config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True
config.General.requestName = "V9_2_W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia"
config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'PAT2CAT_cfg.py'

config.section_("Data")
config.Data.publication  = False
#################################################################
# ALLOWS NON VALID DATASETS
config.Data.allowNonValidInputDataset = True

config.section_("Site")
# Where the output files will be transmitted to
#config.Site.storageSite = 'T2_KR_KNU'
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/CAT/
config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T3_KR_UOS'
config.Data.outLFNDirBase = '/store/group/CAT/V9_2/W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8' 
#config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.blacklist = ['T2_US_Florida']

config.Data.outputDatasetTag='V9_2_RunIIFall17MiniAOD-94X_mc2017_realistic_v11-v1'
config.JobType.pyCfgParams = ['runOnMC=True','useMiniAOD=True','globalTag=94X_mc2017_realistic_v12','runGenTop=False']
config.Data.userInputFiles = lfnList
config.Data.splitting='FileBased'
config.Data.unitsPerJob=1

#config.section_("Debug")
#config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
