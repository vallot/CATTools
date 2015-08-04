import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
]

process.load("CATTools.CatProducer.pseudoTop_cff")

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_Ana",
    )
)

process.outPath = cms.EndPath(process.out)

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("pseudoTop", "partons"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(True),
    status = cms.untracked.vint32( range(1,1000) )
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("pseudoTop", "partons"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(True),
)

process.printList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("pseudoTop", "partons"),
    #src = cms.InputTag("pseudoTop"),
    maxEventsToPrint  = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.partonTop + process.pseudoTop
  #* process.printTree
  #* process.printDecay
  #* process.printList
)

