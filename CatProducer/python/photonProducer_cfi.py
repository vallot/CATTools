import FWCore.ParameterSet.Config as cms

catPhotons = cms.EDProducer("CATPhotonProducer",
    src = cms.InputTag("slimmedPhotons"),
    photonIDs = cms.vstring("cutBasedPhotonID-Spring15-25ns-V1-standalone-loose",
                            "cutBasedPhotonID-Spring15-25ns-V1-standalone-medium",
                            "cutBasedPhotonID-Spring15-25ns-V1-standalone-tight",
                            "mvaPhoID-Spring15-25ns-nonTrig-V2-wp90",
                            ),
    rhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag('catVertex'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    photonIDSources = cms.PSet()

)
