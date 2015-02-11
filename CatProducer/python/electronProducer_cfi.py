import FWCore.ParameterSet.Config as cms

catElectrons = cms.EDProducer("CATElectronProducer",
    src = cms.InputTag("patElectrons"),
    shiftedEnDownSrc = cms.InputTag("shiftedSlimmedElectronsEnDown"),
    shiftedEnUpSrc = cms.InputTag("shiftedSlimmedElectronsEnUp"),
    mcLabel = cms.InputTag("genParticles"),
    vertexLabel = cms.InputTag('offlinePrimaryVertices'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    rhoLabel = cms.InputTag("fixedGridRhoAll", "rho"),
)
