import FWCore.ParameterSet.Config as cms

csvWeights = cms.EDProducer("CSVWeightProducer",
    src = cms.InputTag("catJets"),
    minPt = cms.double(20),
    maxEta = cms.double(2.41),
)
