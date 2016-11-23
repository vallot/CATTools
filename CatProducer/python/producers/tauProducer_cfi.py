import FWCore.ParameterSet.Config as cms

catTaus = cms.EDProducer("CATTauProducer",
    src = cms.InputTag("slimmedTaus"),
)
