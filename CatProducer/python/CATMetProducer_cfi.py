import FWCore.ParameterSet.Config as cms

catMETsSource = "patMETsPFlow"

catMETs = cms.EDProducer("CATMETProducer",
    # input collection
    src = cms.InputTag(catMETsSource),
)

