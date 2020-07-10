import FWCore.ParameterSet.Config as cms

catJets = cms.EDProducer('CATJetProducer',
    src = cms.InputTag('slimmedJets'),
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    btagNames = cms.vstring("pfDeepCSVJetTags:probb","pfDeepCSVJetTags:probbb","pfDeepCSVJetTags:probc", "pfDeepCSVJetTags:probudsg", "pfDeepFlavourJetTags:probb","pfDeepFlavourJetTags:probbb","pfDeepFlavourJetTags:problepb","pfDeepFlavourJetTags:probc","pfDeepFlavourJetTags:probuds","pfDeepFlavourJetTags:probg"),
    payloadName = cms.string('AK4PFchs'),
    jetResFile   = cms.string("CATTools/CatProducer/data/JER/Autumn18_V7b_MC_PtResolution_AK4PFchs.txt"),
    jetResSFFile = cms.string("CATTools/CatProducer/data/JER/Autumn18_V7b_MC_SF_AK4PFchs.txt"),

    setGenParticle = cms.bool(True),
)
