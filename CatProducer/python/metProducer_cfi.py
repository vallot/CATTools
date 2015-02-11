import FWCore.ParameterSet.Config as cms

catMETs = cms.EDProducer("CATMETProducer",
    src = cms.InputTag("patMETsPFlow"),
)
