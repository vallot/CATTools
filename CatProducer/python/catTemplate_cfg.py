import FWCore.ParameterSet.Config as cms

process = cms.Process("CAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
## Standard setup
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

## Options and Output Report
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

## Source
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())

## Max Number of Events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

## total event counter
process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.p = cms.Path(process.nEventsTotal)

## Output Module Configuration (expects a path 'p')
process.catOut = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('catTuple.root'),
    outputCommands = cms.untracked.vstring('drop *')
)

process.schedule = cms.Schedule()    
