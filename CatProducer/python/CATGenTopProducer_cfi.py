import FWCore.ParameterSet.Config as cms

catGenTops = cms.EDProducer("CATGenTopProducer",
    # input collection
    genJetLabel = cms.InputTag("ak5GenJets"),
    mcParticleLabel = cms.InputTag("genParticles"),
)

