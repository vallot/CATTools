import FWCore.ParameterSet.Config as cms

catDstars = cms.EDProducer("CATDStarProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  maxNumPFCand = cms.int32(999),
  d0MassWindow = cms.double(0.05),
  applyCut = cms.bool(True)
)
