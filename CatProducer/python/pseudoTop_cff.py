import FWCore.ParameterSet.Config as cms

partonTop = cms.EDProducer("PartonTopProducer",
    genParticles = cms.InputTag("prunedGenParticles"),
)

from TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi import *
