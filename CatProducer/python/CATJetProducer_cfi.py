import FWCore.ParameterSet.Config as cms

catJets = cms.EDProducer("CATJetProducer",
    # input collection
    src = cms.InputTag("selectedPatJetsPFlow"),
)

