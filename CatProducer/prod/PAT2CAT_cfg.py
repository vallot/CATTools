from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=False # for jpsi candidates
    
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
    process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
    process.out.outputCommands.extend(catEventContentMC)
    if not useMiniAOD:
        process.out.outputCommands.extend(catEventContentAODMC)
if doSecVertex:
    process.out.outputCommands.extend(catEventContentSecVertexs)


####################################################################
#### cmsRun options
####################################################################
process.maxEvents.input = options.maxEvents

# better to have a default file here for test purpose
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       # at CERN
       '/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/MINIAODSIM/MCRUN2_74_V9-v1/00000/2403409D-1225-E511-B64E-0025905A6132.root'
       # AODSIM
       #'/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/MCRUN2_74_V9-v1/00000/54F6E09C-1225-E511-842B-0025905A612E.root'
       # to be added at KISTI
       # '',
    )
)

if options.inputFiles:
    process.source.fileNames = options.inputFiles
#pat input files are removed because it would not work if useMiniAOD is on.    

print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = True 
