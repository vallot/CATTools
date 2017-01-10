import FWCore.ParameterSet.Config as cms
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

catElectrons = cms.EDProducer("CATElectronProducer",
    src = cms.InputTag("slimmedElectrons"),
    unsmaredElectrons = cms.InputTag("slimmedElectrons"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag('catVertex'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    rhoLabel = cms.InputTag("fixedGridRhoAll"),
    minPt = cms.double(5.0),
    maxEta = cms.double(2.5),
    electronIDSources = cms.PSet(),
    electronIDs = cms.vstring(), ## Defined in CatProducer/python/patTools/egmVersionedID_cff.py
)
