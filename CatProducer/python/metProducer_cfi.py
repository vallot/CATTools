import FWCore.ParameterSet.Config as cms

catMETs = cms.EDProducer("CATMETProducer",
    src = cms.InputTag("slimmedMETs"),
)
catMETsPuppi = cms.EDProducer("CATMETProducer",
    src = cms.InputTag("slimmedMETsPuppi"),
)
#There is a CMS rule that we are supposed to use one module per cfi file.  
