import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-0_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_125444/0000/catTuple_1.root',
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.ttljEventSelector_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllGenFilters_cff")
process.load("CATTools.Validation.validation_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

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

process.p = cms.Path(
    process.agen + process.filterPartonTTLJ
  * process.gen + process.rec
  * process.eventsTTLJ
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
