import FWCore.ParameterSet.Config as cms

gen = cms.EDAnalyzer("CATGenValidation",
    weight = cms.InputTag("genWeight"),
    scaleupWeights = cms.InputTag("flatGenWeights:scaleup"),
    scaledownWeights = cms.InputTag("flatGenWeights:scaledown"),
    pdfWeights = cms.InputTag("flatGenWeights:pdf"),
    otherWeights = cms.InputTag("flatGenWeights:other"),
    genParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets"),
)

rec = cms.EDAnalyzer("CATHisAnalysis",
    isMC = cms.untracked.bool(True),
    nPV = cms.InputTag("catVertex:nGoodPV"),
    puWeight = cms.InputTag("pileupWeight"),
    genWeight = cms.InputTag("flatGenWeights"),
    genJets = cms.InputTag("slimmedGenJets"),
    electrons = cms.InputTag("catElectrons"),
    muons = cms.InputTag("catMuons"),
    photons = cms.InputTag("catPhotons"),
    taus = cms.InputTag("catTaus"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
)
  
