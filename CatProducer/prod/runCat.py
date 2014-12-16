from PhysicsTools.PatAlgos.patTemplate_cfg import *
#from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: True default")

options.parseArguments()
runOnMC = options.runOnMC

####################################################################################################
## running PAT
postfix = "PFlow"
jetAlgo="AK5"
from CATTools.CatProducer.catPatSetup_cff import *
catPatConfig(process, runOnMC, postfix, jetAlgo)

process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_DATA_AK5PF'),
            ),
      ), 
      connect = cms.string('sqlite:Winter14_V5_DATA.db')
)

process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
####################################################################################################

from CATTools.CatProducer.catSetup_cff import *
catSetup(process, runOnMC, doSecVertex)

process.maxEvents.input = options.maxEvents

process.source.fileNames = options.inputFiles

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False
