import FWCore.ParameterSet.Config as cms

process = cms.Process("CAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
## Standard setup
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.GlobalTag.globaltag = cms.string( 'START53_V27::All' )

## Options and Output Report
process.options = cms.untracked.PSet(
    #allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

## Source
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())

## Max Number of Events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

## total event counter
process.totaEvents = cms.EDProducer("EventCountProducer")
process.p = cms.Path(process.totaEvents)

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('catTuple.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_cat*_*_*',
        'keep *_goodOfflinePrimaryVertices*_*_*',
        'keep GenEventInfoProduct_*_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep edmMergeableCounter_*_*_*',
        'keep patTriggerPaths_patTrigger*_*_*',
        'keep recoGenParticles_genParticles__SIM',
    )
)
process.outpath = cms.EndPath(process.out)
