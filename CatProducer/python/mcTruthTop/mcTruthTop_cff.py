import FWCore.ParameterSet.Config as cms

partonTop = cms.EDProducer("PartonTopProducer",
    genParticles = cms.InputTag("prunedGenParticles"),
    jetMinPt = cms.double(20),
    jetMaxEta = cms.double(2.5),
    jetConeSize = cms.double(0.4),
)

from TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi import *

from CATTools.CatProducer.mcTruthTop.GenTtbarCategorizer_cfi import *
from CATTools.CatProducer.mcTruthTop.genTopProducer_cfi import *
