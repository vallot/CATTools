import FWCore.ParameterSet.Config as cms

process = cms.Process("h2muAnalyzer")
savename ="h2mu"
datadir = "/xrootd/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150820_215635/0000/"

savename+=".root"

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

import os
for f in os.listdir(datadir):
    process.source.fileNames.append("file:"+datadir)

process.h2mu = cms.EDAnalyzer("h2muAnalyzer",
    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(savename)
)

process.p = cms.Path(process.h2mu)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
