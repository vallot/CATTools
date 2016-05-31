import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/user/jhgoh/CATTools/sync/v7-6-3/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllEventSelector_cfi")
process.load("CATTools.CatAnalyzer.ttll.ttllGenFilters_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllAnalyzers_cff")
process.load("CATTools.CatAnalyzer.ttll.ntuple_cff")
process.ttll.isTopMC = True
process.ttbbll.isTopMC = True

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
process.agen = cms.EDAnalyzer("CATGenTopAnalysis",
    weightIndex = cms.int32(-1),
    weight = cms.InputTag("genWeight"),
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
)

process.p = cms.Path(
    process.agen + process.filterPartonTTLL
  * process.eventsTTLL * process.ttll + process.ttbbll
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
