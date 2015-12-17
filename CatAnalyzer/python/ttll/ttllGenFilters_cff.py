import FWCore.ParameterSet.Config as cms

filterPartonTTLL = cms.EDFilter("TTLLGenCategoryFilter",
    invert = cms.bool(False),
    inputType = cms.string("PartonTop"),
    vetoTau = cms.bool(False),

    src = cms.InputTag("partonTop"),
    addJetType = cms.string("*"),
)

