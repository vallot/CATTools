import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    #'/store/user/jhgoh/CATTools/sync/v7-6-1/DoubleEG_Run2015D-16Dec2015-v2.root',
    #'/store/user/jhgoh/CATTools/sync/v7-6-1/DoubleMuon_Run2015D-16Dec2015-v1.root',
    '/store/user/jhgoh/CATTools/sync/v7-6-1/MuonEG_Run2015D-16Dec2015-v1.root',
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllEventSelector_cfi")
process.load("CATTools.CatAnalyzer.ttll.ttllAnalyzers_cff")
process.load("CATTools.CatAnalyzer.ttll.ntuple_cff")
process.eventsTTLL.isMC = False
if hasattr(process.ntuple.float, "weight"):
    delattr(process.ntuple.float, "weight")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.filterLumi
  * process.eventsTTLL * process.ttll + process.ttbbll
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
