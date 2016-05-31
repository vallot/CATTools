import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
  '/store/group/CAT/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-2_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150923_192710/0000/catTuple_1.root'
]

process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
process.z = cms.EDAnalyzer("CATGenLeptonAnalysis",
    src = cms.InputTag("prunedGenParticles"),
    weight = cms.InputTag("flatGenWeight"),
)

process.p = cms.Path(process.z)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

