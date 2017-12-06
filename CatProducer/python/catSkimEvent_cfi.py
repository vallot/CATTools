import FWCore.ParameterSet.Config as cms

catSkimEvent = cms.EDFilter("CATSkimEventFilter",
    electrons = cms.InputTag("catElectrons"),
    muons = cms.InputTag("catMuons"),
    jets = cms.InputTag("catJets"),
    electronIdNames = cms.vstring(),

    minNLeptons = cms.uint32(0),
    minLeptonPt = cms.double(20),
    maxLeptonAbseta = cms.double(2.4),

    minNJets = cms.uint32(0),
    minJetPt = cms.double(30),
    maxJetAbseta = cms.double(2.5),
)
