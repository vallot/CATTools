from CATTools.CatProducer.catTemplate_cfg import *
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")
options.register('runGenTop', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runGenTop: 1  default")
options.register('runParticleTop', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runParticleTop: 0  default")
options.register('isSignal', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isSignal: 1 default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD
globalTag = options.globalTag
if runOnMC: runGenTop = options.runGenTop
else: runGenTop = False
runParticleTop = False
if runOnMC and options.runParticleTop: runParticleTop = True
isMCSignal = (runOnMC and options.isSignal == True)

####################################################################
#### setting up global tag
####################################################################
#from Configuration.AlCa.autoCond_condDBv2 import autoCond
#if runOnMC: process.GlobalTag.globaltag = autoCond['run2_mc']
#else: process.GlobalTag.globaltag = autoCond['run2_data']
if not globalTag:
    if runOnMC: from CATTools.CatProducer.catDefinitions_cfi import globalTag_mc as globalTag
    else: from CATTools.CatProducer.catDefinitions_cfi import globalTag_rd as globalTag
process.GlobalTag.globaltag = globalTag
print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag    
####################################################################
#### cat tools output
####################################################################
from CATTools.CatProducer.catCandidates_cff import *
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process = addCatCommonObjects(process)
if runOnMC: process = addCatCommonMCObjects(process)
if runGenTop     : process = addCatGenTopObjects(process)
if runParticleTop: process = addCatParticleTopObjects(process)
if isMCSignal:
    process.genWeight.keepFirstOnly = False
    process.catOut.outputCommands.extend(catEventContentMCSignal)

from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeOutput
miniAOD_customizeOutput(process.catOut)
    
####################################################################
#### setting up cat tools
####################################################################
from CATTools.CatProducer.catTools_cff import *
catTool(process, runOnMC, useMiniAOD)
#### add electron ID
addEgmID(process, runOnMC)
####################################################################
#### setting up pat tools - miniAOD step or correcting miniAOD
####################################################################
from CATTools.CatProducer.patTools.patTools_cff import *
patTool(process, runOnMC, useMiniAOD)
#### Finish Paths and Tasks
process.nEventsFiltered = cms.EDProducer("EventCountProducer")
process.p += process.nEventsFiltered
####################################################################
#### cmsRun options
####################################################################
process.maxEvents.input = options.maxEvents

# Default file here for test purpose
if not options.inputFiles:
    if useMiniAOD:
        process.source.fileNames = cms.untracked.vstring('file:/pnfs/knu.ac.kr/data/cms/store/user/jipark/public/cat10x_samples/TTToSemiLeptonic_062A981D-4A57-664A-A583-E803A658594B.root')

else:
    process.source.fileNames = options.inputFiles

#pat input files are removed because it would not work if useMiniAOD is on.    
## to suppress the long output at the end of the job

if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options.wantSummary = True
#process.MessageLogger.cerr.threshold = 'ERROR'
#process.MessageLogger.suppressWarning = cms.untracked.vstring(["JetPtMismatchAtLowPt", "NullTransverseMomentum"])

## for debugging
#process.source.skipEvents = cms.untracked.uint32(3000)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
#print "process.catOut.outputCommands", process.catOut.outputCommands
