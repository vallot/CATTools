import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
if os.path.exists("TTJets_MSDecays_central.txt"):
    mylist = FileUtils.loadListFromFile("TTJets_MSDecays_central.txt")
    process.source.fileNames = mylist
else:
    process.source.fileNames = [
        '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/FEA8DB51-F343-E311-A964-848F69FD2484.root',
    ]

process.pseudoTop = cms.EDProducer("PseudoTopProducer",
    finalStates = cms.InputTag("genParticles"),
    genParticles = cms.InputTag("genParticles"),
    jetConeSize = cms.double(0.5),
    jetMaxEta = cms.double(2.4),
    jetMinPt = cms.double(20),
    leptonConeSize = cms.double(0.1),
    leptonMaxEta = cms.double(2.4),
    leptonMinPt = cms.double(20),
    tMass = cms.double(172.5),
    wMass = cms.double(80.4)
)

process.partonTop = cms.EDProducer("PartonTopProducer",
    genParticles = cms.InputTag("genParticles"),
    jetConeSize = cms.double(0.5),
    jetMaxEta = cms.double(2.4),
    jetMinPt = cms.double(20)
)

process.ttbar = cms.EDAnalyzer("CATGenTopAnalysis",
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
    weight = cms.InputTag("genWeight", "genWeight"),
    weightIndex = cms.uint32(0),
)

process.ttbarNoTau = process.ttbar.clone(filterTaus = cms.bool(True))

process.p = cms.Path(
    process.ttbar + process.ttbarNoTau
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

