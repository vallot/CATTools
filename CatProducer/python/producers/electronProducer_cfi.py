import FWCore.ParameterSet.Config as cms
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

catElectrons = cms.EDProducer("CATElectronProducer",
    src = cms.InputTag("slimmedElectrons"),
    unsmaredElectrons = cms.InputTag("slimmedElectrons"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    vertexLabel = cms.InputTag('catVertex'),
    pfSrc  = cms.InputTag("packedPFCandidates"),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    rhoLabel = cms.InputTag("fixedGridRhoAll"),
    electronIDSources = cms.PSet(),
    electronIDs = cms.vstring(), ## Defined in CatProducer/python/patTools/egmVersionedID_cff.py
    trigObjs = cms.InputTag("slimmedPatTrigger"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
)
