import FWCore.ParameterSet.Config as cms

filterParton = cms.EDFilter("TTGenCategoryFilter",
    invert = cms.bool(False), # set True to invert the selection
    inputType = cms.string("PartonTop"), # category by partonTop
    vetoTau = cms.bool(False), # set True to remove tau->e/mu decay
    nLepton = cms.int32(2), # number of leptons in the final state, negative value to ignore the filter

    src = cms.InputTag("partonTop"),
)

filterGenTop = cms.EDFilter("TTGenCategoryFilter",
    invert = cms.bool(False), # set True to invert the selection
    inputType = cms.string("GenTop"), # category by GenTop object
    vetoTau = cms.bool(False), # set True to remove tau->e/mu decay
    nLepton = cms.int32(2), # number of leptons in the final state, negative value to ignore the filter

    src = cms.InputTag("catGenTops"),
    addJetChannel = cms.string("TTBB"),
)

