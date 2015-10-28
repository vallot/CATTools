import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("h2muAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/xrootd/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-4_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_143012/0000/catTuple_85.root')
    fileNames = cms.untracked.vstring('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/DoubleMuon/v7-4-4_Run2015C_25ns-05Oct2015-v1/151023_165157/0000/catTuple_10.root')

)
        
process.h2mu = cms.EDAnalyzer("h2muAnalyzer",
    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerObjects = cms.InputTag("catTrigger"),
    #triggerObjects = cms.InputTag("selectedPatTrigger"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("h2mu.root")
)

process.p = cms.Path(process.h2mu)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
