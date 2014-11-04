import FWCore.ParameterSet.Config as cms

catGenJetsSource = "ak5GenJets"
catMCsource = "genParticles"

catGenTops = cms.EDProducer("CATGenTopProducer",
    # input collection
    genJetLabel = cms.InputTag(catGenJetsSource),
    mcParticleLabel = cms.InputTag(catMCsource),
)

