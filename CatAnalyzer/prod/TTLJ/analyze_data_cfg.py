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
process.load("CATTools.Validation.ttljEventSelector_cff")
process.load("CATTools.Validation.validation_cff")

process.eventsTTLJ.isMC = False
process.rec.isMC = False
if hasattr(process.eventsTTLJ, "genWeight"): delattr(process.eventsTTLJ, "genWeight")
if hasattr(process, "flatGenWeights"): delattr(process, "flatGenWeights")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.eventsTTLJ.skipHistograms = True
process.eventsTTLJ.applyFilterAt = 1 ## save events from step 1 one lepton

process.load("CATTools.CatAnalyzer.csvWeights_cfi")
process.filterRECO = process.filterRECOMC.clone()
delattr(process, 'filterRECOMC')

from CATTools.CatAnalyzer.analyzers.ntuple_cff import *
process = ntupler_load(process, "eventsTTLJ")

process.pTTLJ = cms.Path(
    process.filterLumi * process.removeLumisWithL1TCert
#  * process.rec
  * process.eventsTTLJ
  * process.ntuple
)

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)

