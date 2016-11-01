import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    #'/store/user/jhgoh/CATTools/sync/v7-6-3/DoubleEG_Run2015D-16Dec2015-v2.root',
    #'/store/user/jhgoh/CATTools/sync/v7-6-3/DoubleMuon_Run2015D-16Dec2015-v1.root',
    '/store/user/jhgoh/CATTools/sync/v7-6-5/MuonEG_Run2015D-16Dec2015-v1.root',
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.ttllEventSelector_cff")
process.load("CATTools.Validation.validation_cff")
process.eventsTTLL.isMC = False
process.rec.isMC = False
if hasattr(process.eventsTTLL, "genWeight"): delattr(process.eventsTTLL, "genWeight")
if hasattr(process, "flatGenWeights"): delattr(process, "flatGenWeights")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.filterLumi
  * process.rec
  * process.eventsTTLL
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
