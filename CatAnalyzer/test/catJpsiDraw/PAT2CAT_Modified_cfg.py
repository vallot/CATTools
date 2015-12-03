from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")
options.register('runGenTop', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runGenTop: 1  default")
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
    process.load("CATTools.CatProducer.genWeight_cff")
    process.load("CATTools.CatProducer.pileupWeight_cff")
    process.catOut.outputCommands.extend(catEventContentMC)
    
if runGenTop:
    from CATTools.CatProducer.catGenHFHadronMatching_cff import *
    genHFTool(process, useMiniAOD)
    process.load("CATTools.CatProducer.mcTruthTop.mcTruthTop_cff")
    process.catOut.outputCommands.extend(catEventContentTOPMC)
    if not useMiniAOD:
        process.catOut.outputCommands.extend(['keep *_catGenTops_*_*',])
        process.catOut.outputCommands.extend(['keep *_muon*_*_*',])
        process.catOut.outputCommands.extend(['keep *_electron*_*_*',])
        process.catOut.outputCommands.extend(['keep *_generalTracks*_*_*',])
        process.catOut.outputCommands.extend(['keep *_*pat*_*_*',])
            
if doSecVertex:
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.catOut.outputCommands.extend(catEventContentSecVertexs)

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
            process.source.fileNames = ['/store/relval/CMSSW_7_4_15/RelValTTbar_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/0253820F-4772-E511-ADD3-002618943856.root']
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
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False

process.catV2Match = cms.EDProducer("MCMatcher",     # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("catSecVertexsV2"),        # RECO objects to match
    matched = cms.InputTag("genParticles"), # mc-truth particle collection
    mcPdgId     = cms.vint32(443),           # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(3),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.15),            # Minimum deltaR for the match
    maxDPtRel = cms.double(0.05),            # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True) # False = just match input in order; True = pick lowest deltaR pair first
)
process.catV3Match = cms.EDProducer("MCMatcher",     # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("catSecVertexsV3"),        # RECO objects to match
    matched = cms.InputTag("genParticles"), # mc-truth particle collection
    mcPdgId     = cms.vint32(443),           # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(3),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.15),            # Minimum deltaR for the match
    maxDPtRel = cms.double(0.05),            # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True) # False = just match input in order; True = pick lowest deltaR pair first
)
process.catV4Match = cms.EDProducer("MCMatcher",     # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("catSecVertexsV4"),        # RECO objects to match
    matched = cms.InputTag("genParticles"), # mc-truth particle collection
    mcPdgId     = cms.vint32(443),           # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(3),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.15),            # Minimum deltaR for the match
    maxDPtRel = cms.double(0.05),            # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True) # False = just match input in order; True = pick lowest deltaR pair first
)

## for debugging
#process.options.wantSummary = True
#process.source.skipEvents = cms.untracked.uint32(3000)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
#print "process.catOut.outputCommands", process.catOut.outputCommands

if doSecVertex :
  if not useMiniAOD :
    process.catOut.outputCommands.extend(['keep *_catSecVertexsV2_*_*',     ])
    process.catOut.outputCommands.extend(['keep *_catV2Match_*_*',     ])
  process.catOut.outputCommands.extend(['keep *_catSecVertexsV3_*_*',     ])
  process.catOut.outputCommands.extend(['keep *_catSecVertexsV4_*_*',     ])
  process.catOut.outputCommands.extend(['keep *_catV3Match_*_*',     ])
  process.catOut.outputCommands.extend(['keep *_catV4Match_*_*',     ])
  process.catOut.outputCommands.extend(['keep *_genParticles_*_*',     ])

## for Jpsi to MuMu samples
process.catOut.outputCommands.remove('keep *_genWeight_*_*')
#process.options.skipEvent=
