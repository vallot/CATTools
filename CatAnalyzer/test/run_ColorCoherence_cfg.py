import FWCore.ParameterSet.Config as cms
process = cms.Process("ColorCoherenceAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-4_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_143547/0000/catTuple_96.root')

process.load("CATTools.CatProducer.pileupWeight_cff")
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.pileupWeight.weightingMethod = "RedoWeight"
process.pileupWeight.pileupRD = pileupWeightMap["Run2015_25nsV1"]
process.pileupWeight.pileupUp = pileupWeightMap["Run2015Up_25nsV1"]
process.pileupWeight.pileupDn = pileupWeightMap["Run2015Dn_25nsV1"]

print process.source.fileNames
process.cc = cms.EDAnalyzer("ColorCoherenceAnalyzer",
    vtx = cms.InputTag("catVertex", "nGoodPV"),
    pileupWeight = cms.InputTag("pileupWeight"),    
    pileupWeight_up = cms.InputTag("pileupWeight","up"),    
    pileupWeight_dn = cms.InputTag("pileupWeight","dn"),    
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("qcd.root")
)

process.p = cms.Path(process.cc)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
