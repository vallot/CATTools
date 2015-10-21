import FWCore.ParameterSet.Config as cms

catMETs = cms.EDProducer("CATMETProducer",
    src = cms.InputTag("slimmedMETs"),
    setUnclusteredEn = cms.bool(True),
)
catMETsPuppi = cms.EDProducer("CATMETProducer",
    src = cms.InputTag("slimmedMETsPuppi"),
    setUnclusteredEn = cms.bool(True),
)
#There is a CMS rule that we are supposed to use one module per cfi file.  
