import FWCore.ParameterSet.Config as cms

catPhotons = cms.EDProducer("CATPhotonProducer",
    # input collection
    src = cms.InputTag("selectedPatPhotons"),
)

