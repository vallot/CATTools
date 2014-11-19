from PhysicsTools.PatAlgos.patTemplate_cfg import *
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
postfix = "PFlow"
jetAlgo="AK5"
doSecVertex=True # for jpsi candidates

from CATTools.CatProducer.catPatSetup_cff import *
from CATTools.CatProducer.catSetup_cff import *
catPatConfig(process, runOnMC, postfix, jetAlgo)
catSetup(process, runOnMC, doSecVertex)

process.maxEvents.input = 100
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source.fileNames = cms.untracked.vstring(
"/store/mc/HC/GenericTTbar/GEN-SIM-RECO/CMSSW_5_3_1_START53_V5-v1/0010/FC85224E-EAAD-E111-AB01-0025901D629C.root")
