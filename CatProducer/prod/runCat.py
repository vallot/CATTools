from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: True default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
globalTag = options.globalTag


from Configuration.AlCa.autoCond import autoCond
if runOnMC:
    process.GlobalTag.globaltag = autoCond['startup']
else:
    process.GlobalTag.globaltag = autoCond['com10']
if globalTag:
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')

print "using globaltag", process.GlobalTag.globaltag

####################################################################################################
## running PAT
postfix = "PFlow"
jetAlgo="AK5"
from CATTools.CatProducer.catPatSetup_cff import *
catPatConfig(process, runOnMC, postfix, jetAlgo)

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
