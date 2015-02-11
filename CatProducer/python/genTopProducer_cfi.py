import FWCore.ParameterSet.Config as cms

catGenTops = cms.EDProducer("CATGenTopProducer",
    # input collection
    genJetLabel = cms.InputTag("ak4GenJets"),
    mcParticleLabel = cms.InputTag("genParticles"),
)

