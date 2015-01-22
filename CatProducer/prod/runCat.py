from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD
globalTag = options.globalTag

print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD

from Configuration.AlCa.GlobalTag import GlobalTag
if useMiniAOD:
    if runOnMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'PLS170_V7AN2::All', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_70_V2_AN1::All', '')

####################################################################################################
## from miniAOD/patTuple_mini.py to run miniAOD maker when starting from AOD
if not useMiniAOD:
    if not globalTag:
        print "ERROR!!!! Need correct globalTag to run on AOD"
    process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')
    process.load('Configuration.StandardSequences.PAT_cff')
    process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import *
    if runOnMC:
        miniAOD_customizeAllMC(process)
    else :
        miniAOD_customizeAllData(process)
####################################################################################################
## setting up catTools
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

from CATTools.CatProducer.catSetup_cff import *
catSetup(process, runOnMC, doSecVertex)

process.maxEvents.input = options.maxEvents
process.source.fileNames = options.inputFiles
## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False
