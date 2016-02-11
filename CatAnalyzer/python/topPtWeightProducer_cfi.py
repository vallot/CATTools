import FWCore.ParameterSet.Config as cms

topPtWeight = cms.EDProducer("TopPtWeightProducer",
    src = cms.InputTag("partonTop"),
)
