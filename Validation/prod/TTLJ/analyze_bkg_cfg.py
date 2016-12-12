import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-3_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/161205_142228/0000/catTuple_1.root',
]

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.Validation.ttljEventSelector_cff")
process.load("CATTools.Validation.validation_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.gen + process.rec
  * process.eventsTTLJ
)

process.load("CATTools.CatProducer.pileupWeight_cff")
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.pileupWeight.weightingMethod = "RedoWeight"
process.pileupWeight.pileupMC = pileupWeightMap["2016_25ns_SpringMC"]
process.pileupWeight.pileupRD = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON"]
process.pileupWeight.pileupUp = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Up"]
process.pileupWeight.pileupDn = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Dn"]
process.eventsTTLJ.vertex.pileupWeight = "pileupWeight::CATeX"

process.eventsTTLJ.filters.filterRECO = "filterRECOMC"
process.eventsTTLJ.filters.ignoreTrig = True

## Customise with cmd arguments
import sys
if len(sys.argv) > 2:
    for l in sys.argv[2:]: exec('process.'+l)
