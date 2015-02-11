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
        'file:../../../catTuple_790.root'
    )
)

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('CATHisAnalysis.root')
)

process.hist = cms.EDAnalyzer("CATHisAnalysis",
    # input collection
    electrons = cms.InputTag("catElectrons"),
    genjets = cms.InputTag("catGenJets"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    muons = cms.InputTag("catMuons"),
    photons = cms.InputTag("catPhotons"),
    #secvertexs = cms.InputTag("catSecVertexs"),
    taus = cms.InputTag("catTaus"),
)


process.p = cms.Path(
            process.hist
)



