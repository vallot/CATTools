import FWCore.ParameterSet.Config as cms

catMETs = cms.EDProducer("CATMETProducer",
    # input collection
    src = cms.InputTag("patMETs"),
)

