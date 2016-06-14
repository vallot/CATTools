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

process.gen = cms.EDAnalyzer("CATGenValidation",
    weight = cms.InputTag("genWeight"),
    scaleupWeights = cms.InputTag("flatGenWeights:scaleup"),
    scaledownWeights = cms.InputTag("flatGenWeights:scaledown"),
    pdfWeights = cms.InputTag("flatGenWeights:pdf"),
    otherWeights = cms.InputTag("flatGenWeights:other"),
    genParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets"),
)

process.p = cms.Path(
    process.gen
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

