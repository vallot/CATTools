import FWCore.ParameterSet.Config as cms

ttll = cms.EDAnalyzer("TTLLAnalyzer",
    reco = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    isTopMC = cms.bool(False),
    bTagName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
)

ttbbll = cms.EDAnalyzer("TTBBLLAnalyzer",
    src = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
)

