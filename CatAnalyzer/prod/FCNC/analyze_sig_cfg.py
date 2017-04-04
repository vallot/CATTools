import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from CATTools.Validation.commonTestInput_cff import commonTestCATTuples
process.source.fileNames = commonTestCATTuples["sig"]
process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.topFCNCEventSelector_cff")
process.load("CATTools.CatAnalyzer.ttll.ttllGenFilters_cff")
process.load("CATTools.Validation.validation_cff")
process.filterParton.nLepton = 1

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
process.load("CATTools.CatProducer.mcTruthTop.partonTop_cfi")
process.agen = cms.EDAnalyzer("CATGenTopAnalysis",
    weightIndex = cms.int32(-1),
    weight = cms.InputTag("flatGenWeights"),
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
)

process.filterTrigEL.triggersToMatch = ["HLT_Ele25_eta2p1_WPTight_Gsf_v", "HLT_Ele27_WPTight_Gsf_v"]
process.eventsFCNC.filters.filterRECO = "filterRECOMC"
process.el = process.eventsFCNC.clone(channel = cms.string("electron"))
process.mu = process.eventsFCNC.clone(channel = cms.string("muon"))
delattr(process, 'eventsFCNC')
#process.el.electron.applyAntiIso = True
#process.mu.muon.applyAntiIso = True

process.load("CATTools.CatAnalyzer.topPtWeightProducer_cfi")

from CATTools.CatAnalyzer.analyzers.ntuple_cff import *
process = ntupler_load(process, "el", "ntupleEL")
process = ntupler_load(process, "mu", "ntupleMU")
process = ntupler_addVarsGen(process, "el", "ntupleEL")
process = ntupler_addVarsGen(process, "mu", "ntupleMU")
process = ntupler_addVarsGenTop(process, "ntupleEL")
process = ntupler_addVarsGenTop(process, "ntupleMU")

process.p_el = cms.Path(
    process.agen + process.filterParton
  * process.gen + process.rec
  * process.el * process.ntupleEL
)

process.p_mu = cms.Path(
    process.agen + process.filterParton
  * process.gen + process.rec
  * process.mu * process.ntupleMU
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
