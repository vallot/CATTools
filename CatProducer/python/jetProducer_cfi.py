import FWCore.ParameterSet.Config as cms

catJets = cms.EDProducer('CATJetProducer',
    src = cms.InputTag('slimmedJets'),
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    btagNames = cms.vstring('pfCombinedInclusiveSecondaryVertexV2BJetTags','pfCombinedMVAV2BJetTags',"inclusiveCandidateSecondaryVerticesCvsL","pfCombinedCvsLJetTags","pfCombinedCvsBJetTags"),
    payloadName = cms.string('AK4PFchs'),
    jetResFile = cms.string("CATTools/CatProducer/data/JER/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt"),
    jetResSFFile = cms.string("CATTools/CatProducer/data/JER/Fall15_25nsV2_MC_SF_AK4PFchs.txt"),
    setGenParticle = cms.bool(True),
)

catJetsPuppi = catJets.clone(
    src = cms.InputTag('slimmedJetsPuppi'),
    btagNames = cms.vstring('pfCombinedInclusiveSecondaryVertexV2BJetTagsAK4PFPuppi','pfCombinedMVAV2BJetTagsAK4PFPuppi'),
    payloadName = cms.string('AK4PFPuppi'),
    jetResFile = cms.string("CATTools/CatProducer/data/JER/Fall15_25nsV2_MC_PtResolution_AK4PFPuppi.txt"),
    jetResSFFile = cms.string("CATTools/CatProducer/data/JER/Fall15_25nsV2_MC_SF_AK4PFPuppi.txt"),
)
#There is a CMS rule that we are supposed to use one module per cfi file.  
