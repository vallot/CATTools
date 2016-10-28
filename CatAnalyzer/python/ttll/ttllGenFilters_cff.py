import FWCore.ParameterSet.Config as cms

filterPartonTTLL = cms.EDFilter("TTGenCategoryFilter",
    invert = cms.bool(False),
    inputType = cms.string("PartonTop"),
    vetoTau = cms.bool(False),
    nLepton = cms.int32(2),

    src = cms.InputTag("partonTop"),
    addJetType = cms.string("*"),
)

filterPartonTTLJ = filterPartonTTLL.clone(
    nLepton = cms.int32(1),
)

