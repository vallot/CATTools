import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from CATTools.Validation.commonTestInput_cff import *
#process.source.fileNames = commonTestCATTuples['sig']
process.source.fileNames = [
    #'/store/group/CAT/ST_FCNC-TH_Tleptonic_HTobb_eta_hut-MadGraph5-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_115237/0000/catTuple_1.root',
    '/store/group/CAT/TT_FCNC-aTtoHJ_Tleptonic_HTobb_eta_hct-MadGraph5-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114948/0000/catTuple_1.root',
]

process.load("CATTools.CatProducer.mcTruthTop.partonTop_cfi")

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_Ana",
    )
)

process.outPath = cms.EndPath(process.out)

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("partonTop", ""),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(True),
    status = cms.untracked.vint32( range(1,1000) )
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("partonTop", ""),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(True),
)

process.printList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("partonTop", ""),
    maxEventsToPrint  = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.partonTop
  #* process.printTree
  #* process.printDecay
  * process.printList
)

