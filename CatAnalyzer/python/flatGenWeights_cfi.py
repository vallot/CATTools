import FWCore.ParameterSet.Config as cms

flatGenWeights = cms.EDProducer("GenWeightsToFlatWeights",
    src = cms.InputTag("genWeight"),
    saveOthers = cms.bool(False),
    keepFirstOnly = cms.bool(True),
)
