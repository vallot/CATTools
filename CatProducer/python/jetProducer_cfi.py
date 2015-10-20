import FWCore.ParameterSet.Config as cms

catJets = cms.EDProducer("CATJetProducer",
    src = cms.InputTag("slimmedJets"),
    btagNames = cms.vstring("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    payloadName = cms.string("AK4PFchs"),
    setGenParticle = cms.bool(True),
)

catJetsPuppi = cms.EDProducer("CATJetProducer",
    src = cms.InputTag("slimmedJetsPuppi"),
    btagNames = cms.vstring("pfCombinedInclusiveSecondaryVertexV2BJetTagsAK4PFPuppi"),
    payloadName = cms.string("AK4PFPuppi"),
    setGenParticle = cms.bool(True),
)
#There is a CMS rule that we are supposed to use one module per cfi file.  
