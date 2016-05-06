import FWCore.ParameterSet.Config as cms

flatGenWeights = cms.EDProducer("GenWeightsFromGenWeight",
    src = cms.InputTag("genWeight")
)
