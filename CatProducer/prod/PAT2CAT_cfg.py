from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
doDstar=True      # for Dstar meson.
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")
options.register('runGenTop', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runGenTop: 1  default")
options.register('runOnRelVal', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnRelVal: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD
globalTag = options.globalTag
if runOnMC: runGenTop = options.runGenTop
else: runGenTop = False

####################################################################
#### setting up global tag
####################################################################
from Configuration.AlCa.autoCond_condDBv2 import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc']
if not runOnMC:
    process.GlobalTag.globaltag = autoCond['run2_data']
if globalTag:
    process.GlobalTag.globaltag = globalTag
print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag    
####################################################################
#### cat tools output
####################################################################
process.load("CATTools.CatProducer.catCandidates_cff")    
from CATTools.CatProducer.catEventContent_cff import *
process.catOut.outputCommands = catEventContent

if runOnMC:
    process.load("CATTools.CatProducer.pileupWeight_cff")
    process.load("CATTools.CatProducer.genWeight_cff")
    process.catOut.outputCommands.extend(catEventContentMC)
else: 
    process.catOut.outputCommands.extend(catEventContentRD)
    
if runGenTop:
    process.load("CATTools.CatProducer.mcTruthTop.mcTruthTop_cff")
    process.catOut.outputCommands.extend(catEventContentTOPMC)
    # for GenTtbarCategories
    from CATTools.CatProducer.Tools.tools import *
    genHFTool(process, useMiniAOD)
    process.catOut.outputCommands.extend(['keep *_catGenTops_*_*',])
            
if doSecVertex:
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.catOut.outputCommands.extend(catEventContentSecVertexs)

if doDstar :
    process.catOut.outputCommands.extend(['keep *_catDstars_*_*',])


from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeOutput
miniAOD_customizeOutput(process.catOut)
    
process.catOutpath = cms.EndPath(process.catOut)    
process.schedule.append(process.catOutpath)
####################################################################
#### setting up cat tools
####################################################################
from CATTools.CatProducer.catTools_cff import *
catTool(process, runOnMC, useMiniAOD)
####################################################################
#### setting up pat tools - miniAOD step or correcting miniAOD
####################################################################
from CATTools.CatProducer.patTools_cff import *
patTool(process, runOnMC, useMiniAOD)
####################################################################
#### cmsRun options
####################################################################
process.maxEvents.input = options.maxEvents

# Default file here for test purpose
if not options.inputFiles:
    if useMiniAOD:
        if runGenTop:
            process.source.fileNames = [
                                      '/store/relval/CMSSW_7_4_15/RelValTTbar_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/0253820F-4772-E511-ADD3-002618943856.root',
                                      '/store/relval/CMSSW_7_4_15/RelValTTbar_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/3848030E-4772-E511-AFCA-0026189438CC.root']
        else:
            process.source.fileNames = ['/store/relval/CMSSW_7_4_15/RelValZMM_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/10FF6E32-3C72-E511-87AD-0025905A60B4.root']
    else:
        if useGenTop:
            process.source.fileNames = ['/store/relval/CMSSW_7_4_12/RelValTTbar_13/GEN-SIM-RECO/PU25ns_74X_mcRun2_asymptotic_v2_v2-v1/00000/006F3660-4B5E-E511-B8FD-0025905B8596.root']
        else:
            process.source.fileNames = ['/store/relval/CMSSW_7_4_15/RelValZMM_13/GEN-SIM-RECO/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/18B13146-3872-E511-A382-00261894394D.root']

if options.inputFiles:
    process.source.fileNames = options.inputFiles
#pat input files are removed because it would not work if useMiniAOD is on.    

## to suppress the long output at the end of the job



if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = True

process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.suppressWarning = cms.untracked.vstring(["JetPtMismatchAtLowPt", "NullTransverseMomentum"])


## for debugging
#process.options.wantSummary = True
#process.source.skipEvents = cms.untracked.uint32(3000)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
#print "process.catOut.outputCommands", process.catOut.outputCommands
