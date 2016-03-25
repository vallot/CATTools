import FWCore.ParameterSet.Config as cms
process = cms.Process("CtagAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['file:catTuple.root']
process.source.fileNames = ['/store/user/geonmo/for2016SpringKPS_v10/TT_TuneCUETP8M1_13TeV-powheg-pythia8/dstar_v11_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160322_070040/0002/catTuple_2413.root']



process.cattree = cms.EDAnalyzer("CATCTagAnalyzer",
    jets = cms.InputTag("catJets"),
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ctag.root"
))

process.p = cms.Path(process.cattree)
#process.MessageLogger.cerr.FwkReport.reportEvery = 50000
