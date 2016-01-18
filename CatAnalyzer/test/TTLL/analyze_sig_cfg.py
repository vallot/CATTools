import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151127_200613/0000/catTuple_1.root'
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllEventSelector_cfi")
process.load("CATTools.CatAnalyzer.ttll.ttllGenFilters_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllAnalyzers_cff")
process.load("CATTools.CatAnalyzer.ttll.ntuple_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.agen = cms.EDAnalyzer("CATGenTopAnalysis",
    weightIndex = cms.uint32(0),
    weight = cms.InputTag("genWeight", "genWeight"),
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
)

process.p = cms.Path(
    process.agen
  + process.filterPartonTTLL * process.eventsTTLL * process.ttbbll
#  * process.ntuple
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
