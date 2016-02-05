import FWCore.ParameterSet.Config as cms

catDstars = cms.EDProducer("CATDStarProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  maxNumPFCand = cms.int32(5),
  applyCut = cms.bool(False)
)
