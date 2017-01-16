import FWCore.ParameterSet.Config as cms

ttLL = cms.EDAnalyzer("TTLLAnalyzer",
    isMC = cms.bool(True),
    isTTbar = cms.bool(False),
    doGenWeightSysts = cms.bool(False),

    nVertex = cms.InputTag("catVertex:nGoodPV"),

    src = cms.InputTag("eventsTTLJ"),

    topPtWeight = cms.InputTag("topPtWeight"),
    genWeight = cms.InputTag("flatGenWeights"),
    pdfWeights = cms.InputTag("flatGenWeights:pdf"),
    scaleWeights = cms.InputTag("flatGenWeights:scale"),
    csvWeight = cms.InputTag("csvWeights"),

    puWeight = cms.InputTag("pileupWeight"),

    partonTop = cms.InputTag("partonTop"),
)

