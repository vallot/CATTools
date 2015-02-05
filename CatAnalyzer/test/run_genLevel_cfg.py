import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
#    'file:/store1/jhgoh/catTuple_790.root',
    'file:/afs/cern.ch/user/c/chanwook/catTuple_290.root',
]

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.partontop = cms.EDAnalyzer("CATGenLevelAnalysis",
    channel = cms.InputTag("partonTop", "channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partons = cms.InputTag("partonTop"),
    leptonMinPt = cms.double(20.),
    leptonMaxEta = cms.double(2.4),
    jetMinPt = cms.double(20.),
    jetMaxEta = cms.double(2.4),
)

process.pseudotop = cms.EDAnalyzer("CATGenLevelAnalysis",
    channel = cms.InputTag("partonTop", "channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partons = cms.InputTag("pseudoTop"),
    leptonMinPt = cms.double(20.),
    leptonMaxEta = cms.double(2.4),
    jetMinPt = cms.double(20.),
    jetMaxEta = cms.double(2.4),
)


process.p = cms.Path(
    process.partontop+process.pseudotop
)


