from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=False # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD
globalTag = options.globalTag

if globalTag:
    process.GlobalTag.globaltag = globalTag

from CATTools.CatProducer.patTools_cff import *
patTool(process, runOnMC, useMiniAOD)

from CATTools.CatProducer.catTools_cff import *
catTool(process, runOnMC, doSecVertex, useMiniAOD)

from CATTools.CatProducer.catEventContent_cff import *
process.out.outputCommands = catEventContent
if runOnMC:
    process.out.outputCommands.extend(catEventContentMC)
    if not useMiniAOD:
        process.out.outputCommands.extend(catEventContentAODMC)
    if doSecVertex:
        process.out.outputCommands.extend(catEventContentSecVertexs)
        
process.maxEvents.input = options.maxEvents
if options.inputFiles:
    process.source.fileNames = options.inputFiles
else :
    from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValTTbarPileUpGENSIMRECO
    process.source.fileNames = filesRelValTTbarPileUpGENSIMRECO
    
print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False
