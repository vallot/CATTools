import FWCore.ParameterSet.Config as cms

catGenJets = cms.EDProducer("CATGenJetProducer",
    # input collection
    src = cms.InputTag("ak5GenJets"),
    pt = cms.double(10),
    eta = cms.double(2.5)
)

