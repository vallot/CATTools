import FWCore.ParameterSet.Config as cms

catPhotons = cms.EDProducer("CATPhotonProducer",
    src = cms.InputTag("slimmedPhotons"),
    unsmearedPhotons = cms.InputTag("slimmedPhotons"),
    rhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag('catVertex'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    photonIDSources = cms.PSet(),
    photonIDs = cms.vstring(),  ## Defined in CatProducer/python/patTools/egmVersionedID_cff.py
)
