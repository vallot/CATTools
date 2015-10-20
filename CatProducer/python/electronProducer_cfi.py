import FWCore.ParameterSet.Config as cms

catElectrons = cms.EDProducer("CATElectronProducer",
    src = cms.InputTag("slimmedElectrons"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag('catVertex'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    rhoLabel = cms.InputTag("fixedGridRhoAll"),
    electronIDSources = cms.PSet(),
    electronIDs = cms.vstring(
        "cutBasedElectronID-Spring15-25ns-V1-standalone-loose",
        "cutBasedElectronID-Spring15-25ns-V1-standalone-medium",
        "cutBasedElectronID-Spring15-25ns-V1-standalone-tight",
        "cutBasedElectronID-Spring15-25ns-V1-standalone-veto",
        "heepElectronID-HEEPV60",
        "mvaEleID-Spring15-25ns-nonTrig-V1-wp80",
        "mvaEleID-Spring15-25ns-nonTrig-V1-wp90"
        )    
)
