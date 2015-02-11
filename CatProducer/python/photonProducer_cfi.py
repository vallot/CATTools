import FWCore.ParameterSet.Config as cms

catPhotons = cms.EDProducer("CATPhotonProducer",
    src = cms.InputTag("selectedPatPhotons"),
)
