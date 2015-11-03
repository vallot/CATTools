import FWCore.ParameterSet.Config as cms
process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("CATTools.CatProducer.pileupWeight_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:catTuple.root'))
process.pileupWeight.weightingMethod = "RedoWeight"

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_pileupWeight*_*_*",
    )
)
process.outPath = cms.EndPath(process.out)
