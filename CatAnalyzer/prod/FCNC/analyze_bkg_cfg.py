import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from CATTools.Validation.commonTestInput_cff import commonTestCATTuples
process.source.fileNames = commonTestCATTuples["bkg"]
process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.topFCNCEventSelector_cff")
process.load("CATTools.Validation.validation_cff")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.filterTrigEL.triggersToMatch = ["HLT_Ele25_eta2p1_WPTight_Gsf_v", "HLT_Ele27_WPTight_Gsf_v"]
process.eventsFCNC.filters.filterRECO = "filterRECOMC"
process.el = process.eventsFCNC.clone(channel = cms.string("electron"))
process.mu = process.eventsFCNC.clone(channel = cms.string("muon"))
delattr(process, 'eventsFCNC')
#process.el.electron.applyAntiIso = True
#process.mu.muon.applyAntiIso = True

from CATTools.CatAnalyzer.analyzers.ntuple_cff import *
process = ntupler_load(process, "el", "ntupleEL")
process = ntupler_load(process, "mu", "ntupleMU")
process = ntupler_addVarsGen(process, "el", "ntupleEL")
process = ntupler_addVarsGen(process, "mu", "ntupleMU")
#process = ntupler_addVarsGenTop(process, "el", "ntupleEL")
#process = ntupler_addVarsGenTop(process, "mu", "ntupleMU")

process.p_el = cms.Path(
    process.gen + process.rec
  * process.el * process.ntupleEL
)

process.p_mu = cms.Path(
    process.gen + process.rec
  * process.mu * process.ntupleMU
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)

process.source.fileNames = ["/store/group/CAT/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_103306/0000/catTuple_1.root"]
