## import skeleton process
import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:catTuple.root'
    )
)

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('CATDStarAnalysis.root')
)

process.defaultIso = cms.EDAnalyzer("CATDStarAnalysis",
    # input collection
    D0Src    = cms.InputTag("catDstars","D0Cand"),
    DstarSrc = cms.InputTag("catDstars","DstarCand"),
)


process.p = cms.Path(
            process.defaultIso
)



