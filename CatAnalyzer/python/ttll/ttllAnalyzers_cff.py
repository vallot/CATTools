import FWCore.ParameterSet.Config as cms
from CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionProducer_cfi import *

ttll = cms.EDAnalyzer("TTLLAnalyzer",
    recoObjects = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    doTree = cms.bool(True),
    isTopMC = cms.bool(False),
    bTagName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    kinfit = cms.InputTag("ttbarDileptonKin"),
)

ttbbll = cms.EDAnalyzer("TTBBLLAnalyzer",
    src = cms.InputTag("eventsTTLL"),
    partonTop = cms.InputTag("partonTop"),
    doTree = cms.bool(True),
    isTopMC = cms.bool(False),
)

