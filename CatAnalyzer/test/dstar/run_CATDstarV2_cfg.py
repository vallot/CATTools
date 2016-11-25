import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonSVAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("CATTools.Validation.ttllEventSelector_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)


## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('isTT',False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isTT: 0  default")
options.parseArguments()

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TTTo2L2Nu_13TeV-powheg/v8-0-0_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/160705_215520/0000/catTuple_1.root']
#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TTTo2L2Nu_13TeV-powheg/v8-0-2_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/161104_204219/0000/catTuple_4.root']
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v8-0-2_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/161104_204003/0000/catTuple_4.root']


#process.MessageLogger.debugModules = cms.untracked.vstring('cattree')
#process.MessageLogger.destinations = cms.untracked.vstring('detailInfo')
#process.MessageLogger.detailInfo = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG'))

from CATTools.CatAnalyzer.leptonSF_cff import *

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllGenFilters_cff")

process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
process.agen = cms.EDAnalyzer("CATGenTopAnalysis",
    weightIndex = cms.int32(-1),
    weight = cms.InputTag("flatGenWeights"),
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
)

process.cattree = cms.EDAnalyzer("CATDstarV2Analyzer",
    recoObjects = cms.InputTag("eventsTTLL"),
    vertices = cms.InputTag("catVertex"),
    d0s    = cms.InputTag("catDstars","D0Cand"),
    dstars = cms.InputTag("catDstars","DstarCand"),
    Jpsis  = cms.InputTag("catDstars","JpsiCand"),
    matchingDeltaR = cms.double(0.15),
    mcLabel = cms.InputTag("prunedGenParticles"),
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("dstar_v2_result.root"
))


#process.eventsTTLL.filters.ignoreTrig = cms.bool(True)
process.eventsTTLL.filters.ignoreTrig = cms.bool(False)
process.eventsTTLL.applyFilterAt = cms.int32(5)


if ( options.isTT ) : 
  print "This is TT Samples. Run agen and filter parto."
  process.p = cms.Path(
      process.agen + process.filterPartonTTLL* process.eventsTTLL * process.cattree
  )
else : 
  process.p = cms.Path( process.eventsTTLL * process.cattree)


if ( process.maxEvents.input <0 or process.maxEvents > 5000) :
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = True




