import FWCore.ParameterSet.Config as cms

catGenTops = cms.EDProducer("CATGenTopProducer",
    # input collection
    genJetLabel = cms.InputTag("ak4GenJets"),
    mcParticleLabel = cms.InputTag("genParticles"),
    genBHadJetIndex = cms.InputTag("matchGenBHadron", "genBHadJetIndex"),
    genBHadFlavour = cms.InputTag("matchGenBHadron", "genBHadFlavour"),
    genCHadJetIndex = cms.InputTag("matchGenCHadron", "genCHadJetIndex"),
    genCHadFlavour = cms.InputTag("matchGenCHadron", "genCHadFlavour"),
)

