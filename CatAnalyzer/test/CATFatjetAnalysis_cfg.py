## import skeleton process
import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:ttbar2.root'
    )
)

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('CATFatjetAnalysis.root')
)

process.defaultIso = cms.EDAnalyzer("CATFatjetAnalysis",
    # input collection
    muonLabel = cms.InputTag("catMuons"),
    electronLabel = cms.InputTag("catElectrons"),
    jetLabel = cms.InputTag("catJets"),
    fatjetLabel = cms.InputTag("catFatJets"),
    metLabel = cms.InputTag("catMETs"),
    Verbosity = cms.untracked.int32(0) # set to 1 (or greater)  for printouts
)

process.MessageLogger.suppressInfo = cms.untracked.vstring('CATFatjetAnalysis')
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.p = cms.Path(
            process.defaultIso
)



