import FWCore.ParameterSet.Config as cms

genJetHadronFlavour = cms.EDProducer("GenJetHadronMatchProducer",
    genJet = cms.InputTag("slimmedGenJets"),
    genParticle = cms.InputTag("prunedGenParticles"),
)
