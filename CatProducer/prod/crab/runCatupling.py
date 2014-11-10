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
process.source.fileNames = cms.untracked.vstring(
'file:/pnfs/user/kraft_data/FEEEC639-4A98-E211-BE1C-002618943919.root',
#'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/001A0DC8-C313-E211-BCCB-00261894397B.root'
#'file:/pnfs/user/qcd/QCD_Pt-15to3000_Tune1_Flat_8TeV_pythia8_AODSIM_PU_S7_START52_V9-v1/02F7EBD9-A09E-E111-AF77-003048C692C0.root'
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/003F2216-14E1-E111-960B-003048C693EE.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/020431EA-28E1-E111-B959-0030487E4ED5.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/02B56FF3-11E1-E111-ABE0-002590494C94.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/02F44798-2DE1-E111-83F3-00266CF327C4.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/088B3B16-38E1-E111-99E6-0030487D5DB1.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/08BEF948-28E1-E111-AAF6-0025901D490C.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/0A2B5181-33E1-E111-8A7C-002590494C74.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/0C3D2D22-43E1-E111-A107-00266CF9C1AC.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/0C608324-36E1-E111-8708-003048F0E5A4.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/14085CA3-3DE1-E111-BB95-00266CF270A8.root'
)
