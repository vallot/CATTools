import FWCore.ParameterSet.Config as cms

ntupleFCNC = cms.EDAnalyzer("FCNCNtupler",
    isMC = cms.bool(True),

    nVertex = cms.InputTag("catVertex:nGoodPV"),

    src = cms.InputTag("event"),

    genWeight = cms.InputTag("flatGenWeights"),
    pdfWeights = cms.InputTag("flatGenWeights:pdf"),
    scaleWeights = cms.InputTag("flatGenWeights:scale"),
    csvWeight = cms.InputTag("csvWeights"),
    csvWeightSyst = cms.InputTag("csvWeights:syst"),

    puWeight = cms.InputTag("pileupWeight"),
)

