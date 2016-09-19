import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/group/CAT/TTbarXSecSynchronization/v8-0-1/TT_TuneCUETP8M1_13TeV-powheg-pythia8__PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext3-v1__1671FA99-240F-E611-BF05-00266CFAE464.root'
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.ttllEventSelector_cfi")
process.load("CATTools.Validation.validation_cff")

eventsTTLL.electron.idName = "cutBasedElectronID-Spring15-25ns-V1-standalone-medium"
eventsTTLL.electron.applyEcalCrackVeto = True
eventsTTLL.jet.bTagName = "pfCombinedInclusiveSecondaryVertexV2BJetTags"
eventsTTLL.jet.bTagWP = "CSVM"
eventsTTLL.jet.skipJER = True
eventsTTLL.filters.ignoreTrig = True

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.gen + process.rec
  * process.eventsTTLL
)

