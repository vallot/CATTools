import FWCore.ParameterSet.Config as cms
process = cms.Process("ColorCoherenceAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
process.source.fileNames.append("file:/pnfs/user/jlee/test/catTuple.root")
print process.source.fileNames
process.cc = cms.EDAnalyzer("ColorCoherenceAnalyzer",
    vtx = cms.InputTag("catVertex", "nPV"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("test_2p76.root")
)

process.p = cms.Path(process.cc)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
