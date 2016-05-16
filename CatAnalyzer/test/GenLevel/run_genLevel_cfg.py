import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/user/jhgoh/CATTools/sync/v7-6-3/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',
]

#process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
#process.load("CATTools.CatProducers.mcTruthTop.partonTop_cfi")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")

process.ttbar = cms.EDAnalyzer("CATGenTopAnalysis",
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
    weight = cms.InputTag("flatGenWeights"),
    weightIndex = cms.int32(-1),
)

process.ttbarNoTau = process.ttbar.clone(filterTaus = cms.bool(True))

process.p = cms.Path(
    process.ttbar + process.ttbarNoTau
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

