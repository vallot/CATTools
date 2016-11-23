import FWCore.ParameterSet.Config as cms

catGenJets = cms.EDProducer("CATGenJetProducer",
    src = cms.InputTag("slimmedGenJets"),
    pt = cms.double(10),
    eta = cms.double(2.5)
)
