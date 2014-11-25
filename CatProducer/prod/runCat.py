from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=False # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD

print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD
####################################################################################################
## from miniAOD/patTuple_mini.py to run miniAOD maker when starting from AOD
if not useMiniAOD:
    process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
    process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
    process.load("PhysicsTools.PatAlgos.slimming.slimming_cff")
    process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")

    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeCommon,miniAOD_customizeMC,miniAOD_customizeData
    miniAOD_customizeCommon(process)
    if runOnMC:
        miniAOD_customizeMC(process)
    else :
        miniAOD_customizeData(process)
####################################################################################################
## setting up catTools
from CATTools.CatProducer.catPatSetup_cff import *
catPatConfig(process, runOnMC)

from CATTools.CatProducer.catSetup_cff import *
catSetup(process, runOnMC, doSecVertex)

process.maxEvents.input = 5000

process.source.fileNames = options.inputFiles

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False
