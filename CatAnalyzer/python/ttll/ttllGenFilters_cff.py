import FWCore.ParameterSet.Config as cms

filterParton = cms.EDFilter("TTGenCategoryFilter",
    invert = cms.bool(False),
    inputType = cms.string("PartonTop"),
    vetoTau = cms.bool(False),
    nLepton = cms.int32(2),

    src = cms.InputTag("partonTop"),
)

filterGenTop = cms.EDFilter("TTGenCategoryFilter",
    invert = cms.bool(False),
    inputType = cms.string("GenTop"),
    vetoTau = cms.bool(False),
    nLepton = cms.int32(2),

    src = cms.InputTag("catGenTops"),
    addJetChannel = cms.string("TTBB"),
)

