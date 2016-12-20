import FWCore.ParameterSet.Config as cms

csvWeights = cms.EDProducer("CSVWeightProducer",
    src = cms.InputTag("catJets"),
    minPt = cms.double(15),
    maxEta = cms.double(2.5),
)
