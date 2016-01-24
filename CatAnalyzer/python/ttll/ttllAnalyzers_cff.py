import FWCore.ParameterSet.Config as cms
from CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionProducer_cfi import *

ttll = cms.EDAnalyzer("TTLLAnalyzer",
    recoObjects = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    isTopMC = cms.bool(False),
    bTagName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    kinfit = cms.InputTag("ttbarDileptonKin"),
)

ttbbll = cms.EDAnalyzer("TTBBLLAnalyzer",
    src = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
    isTopMC = cms.bool(False),
)

