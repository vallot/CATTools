import FWCore.ParameterSet.Config as cms

catElectrons = cms.EDProducer("CATElectronProducer",
    src = cms.InputTag("slimmedElectrons"),
    shiftedEnDownSrc = cms.InputTag("shiftedSlimmedElectronsEnDown"),
    shiftedEnUpSrc = cms.InputTag("shiftedSlimmedElectronsEnUp"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    rhoLabel = cms.InputTag("fixedGridRhoAll", "rho"),
)
