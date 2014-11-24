import FWCore.ParameterSet.Config as cms
process = cms.Process("CAT")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
process.load("Configuration.StandardSequences.MagneticField_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
                                fileName = cms.untracked.string('catTuple.root'),
                                outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning )
                               )
process.outpath = cms.EndPath(process.out)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#from PhysicsTools.PatAlgos.patTemplate_cfg import *
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('runOnMC', True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "runOnMC: 1  default")

import sys
if hasattr(sys, "argv") == True:
    options.parseArguments()
    runOnMC = options.runOnMC

print "runOnMC",runOnMC
useMiniAOD = True
doSecVertex=False # for jpsi candidates
postfix = "PFlow"
jetAlgo="AK5"

from CATTools.CatProducer.catPatSetup_cff import *
catPatConfig(process, runOnMC, postfix, jetAlgo)

from CATTools.CatProducer.catSetup_cff import *
catSetup(process, runOnMC, useMiniAOD, doSecVertex)

process.maxEvents.input = 3000
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/pnfs/user/jlee/DYJetsToLL_M-50_13TeV-madgraph-pythia8/0432E62A-7A6C-E411-87BB-002590DB92A8.root")
)
