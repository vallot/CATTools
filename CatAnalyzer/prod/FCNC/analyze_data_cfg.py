import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from CATTools.Validation.commonTestInput_cff import commonTestCATTuples
process.source.fileNames = commonTestCATTuples["data"]
process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.topFCNCEventSelector_cff")
process.load("CATTools.Validation.validation_cff")
process.eventsFCNC.isMC = False
process.rec.isMC = False
if hasattr(process.eventsFCNC, "genWeight"): delattr(process.eventsFCNC, "genWeight")
if hasattr(process, "flatGenWeights"): delattr(process, "flatGenWeights")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.filterTrigEL.triggersToMatch = ["HLT_Ele25_eta2p1_WPTight_Gsf_v", "HLT_Ele27_WPTight_Gsf_v"]
process.el = process.eventsFCNC.clone(channel = cms.string("electron"))
process.mu = process.eventsFCNC.clone(channel = cms.string("muon"))
delattr(process, 'eventsFCNC')
#process.el.electron.applyAntiIso = True
#process.mu.muon.applyAntiIso = True

from CATTools.CatAnalyzer.analyzers.ntuple_cff import *
process = ntupler_load(process, "el", "ntupleEL")
process = ntupler_load(process, "mu", "ntupleMU")
#process = ntupler_addVarsGen(process, "el", "ntupleEL")
#process = ntupler_addVarsGen(process, "mu", "ntupleMU")
#process = ntupler_addVarsGenTop(process, "el", "ntupleEL")
#process = ntupler_addVarsGenTop(process, "mu", "ntupleMU")

process.p_el = cms.Path(
    process.filterLumi# * process.removeLumisWithL1TCert
  * process.rec
  * process.el * process.ntupleEL
)

process.p_mu = cms.Path(
    process.filterLumi# * process.removeLumisWithL1TCert
  * process.rec
  * process.mu * process.ntupleMU
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
