from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=False # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")
options.register('runGenTop', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runGenTop: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD
globalTag = options.globalTag
if runOnMC: runGenTop = options.runGenTop
else: runGenTop = False

####################################################################
#### setting up global tag
####################################################################
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc_FULL']
if not runOnMC:
    process.GlobalTag.globaltag = autoCond['run2_data']
if globalTag:
    process.GlobalTag.globaltag = globalTag

####################################################################
#### setting up pat tools - miniaod step
####################################################################
from CATTools.CatProducer.patTools_cff import *
patTool(process, runOnMC, useMiniAOD)

####################################################################
#### setting up cat tools
####################################################################
from CATTools.CatProducer.catTools_cff import *
catTool(process, runOnMC, doSecVertex, useMiniAOD)

from CATTools.CatProducer.catEventContent_cff import *
process.out.outputCommands = catEventContent
if runOnMC:
    from CATTools.CatProducer.catGenHFHadronMatching_cff import *
    genHFTool(process, useMiniAOD)
    if runGenTop:
        process.load("CATTools.CatProducer.mcTruthTop.mcTruthTop_cff")
        if not useMiniAOD:
            process.out.outputCommands.extend(catEventContentAODMC)
    process.out.outputCommands.extend(catEventContentMC)
if doSecVertex:
    process.out.outputCommands.extend(catEventContentSecVertexs)


####################################################################
#### cmsRun options
####################################################################
process.maxEvents.input = options.maxEvents

# Default file here for test purpose
if useMiniAOD:
    process.source.fileNames = ['/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/MINIAODSIM/MCRUN2_74_V9-v1/00000/2403409D-1225-E511-B64E-0025905A6132.root']
    ## Hack to run on relval sample
    process.genMetExtractor.metSource = "slimmedMETs::RECO"
else:
    process.source.fileNames = ['/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/MCRUN2_74_V9-v1/00000/54F6E09C-1225-E511-842B-0025905A612E.root']

if options.inputFiles:
    process.source.fileNames = options.inputFiles
#pat input files are removed because it would not work if useMiniAOD is on.    

print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False

## for debugging
#process.source.skipEvents = cms.untracked.uint32(3000)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
