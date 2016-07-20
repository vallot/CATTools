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
    'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-6-5_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160524_085454/0000/catTuple_1.root'
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

