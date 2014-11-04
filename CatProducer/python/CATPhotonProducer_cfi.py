import FWCore.ParameterSet.Config as cms

catPhotonsSource = "selectedPatPhotons"

catPhotons = cms.EDProducer("CATPhotonProducer",
    # input collection
    src = cms.InputTag(catPhotonsSource),
)

