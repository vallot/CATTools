import FWCore.ParameterSet.Config as cms

ttbbll = cms.EDAnalyzer("TTBBLLAnalyzer",
    src = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
)

