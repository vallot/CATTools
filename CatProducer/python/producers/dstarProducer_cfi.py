import FWCore.ParameterSet.Config as cms

catDstars = cms.EDProducer("CATDStarProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  mcLabel  = cms.InputTag("prunedGenParticles"),
  maxNumPFCand = cms.int32(9999),
  maxDeltaR = cms.double(0.2),
  d0MassCut = cms.double(0.5),
  d0MassWindow = cms.double(0.05),
  vertexLabel = cms.InputTag("catVertex"),
  applySoftLeptonCut = cms.bool(True)
)
