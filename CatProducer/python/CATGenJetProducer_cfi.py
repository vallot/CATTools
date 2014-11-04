import FWCore.ParameterSet.Config as cms

catGenJetsSource = "ak5GenJets"

catGenJets = cms.EDProducer("CATGenJetProducer",
    # input collection
    src = cms.InputTag(catGenJetsSource),
    pt = cms.double(10),
    eta = cms.double(2.5)
)

