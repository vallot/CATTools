import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


mylist = FileUtils.loadListFromFile ("~/ttbar_miniaod.list")
readFiles = cms.untracked.vstring( *mylist)
process.source = cms.Source("PoolSource", fileNames = readFiles)
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames = [
#    'file:/afs/cern.ch/user/c/chanwook/catTuple_290.root',
#]

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.ana = cms.EDAnalyzer("CATGenLevelAnalysis",
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partons = cms.InputTag("partonTop"),
    pseudo = cms.InputTag("pseudoTop"),
    leptonMinPt = cms.double(20.),
    leptonMaxEta = cms.double(2.4),
    jetMinPt = cms.double(20.),
    jetMaxEta = cms.double(2.4),
)

process.p = cms.Path(
    process.ana
)
