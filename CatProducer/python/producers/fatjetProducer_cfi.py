import FWCore.ParameterSet.Config as cms

catFatJets = cms.EDProducer('CATFatJetProducer',
    src = cms.InputTag('slimmedJetsAK8'),
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    btagNames = cms.vstring('pfBoostedDoubleSecondaryVertexAK8BJetTags', 'pfCombinedInclusiveSecondaryVertexV2BJetTags',"pfCombinedCvsLJetTags","pfCombinedCvsBJetTags"),
    payloadName = cms.string('AK8PFchs'),
    jetResFile   = cms.string("CATTools/CatProducer/data/JER/Fall17_V3_MC_PtResolution_AK8PFchs.txt"),
    jetResSFFile = cms.string("CATTools/CatProducer/data/JER/Fall17_V3_MC_SF_AK8PFchs.txt"),
    setGenParticle = cms.bool(False),
)

#There is a CMS rule that we are supposed to use one module per cfi file.  
