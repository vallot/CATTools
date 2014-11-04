import FWCore.ParameterSet.Config as cms

catJetsSource = "selectedPatJetsPFlow"

catJets = cms.EDProducer("CATJetProducer",
    # input collection
    src = cms.InputTag(catJetsSource),
)

