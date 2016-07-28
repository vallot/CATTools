import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/MuonEG_Run2015D-16Dec2015-v1.root',]
#process.source.fileNames = ['file:/xrootd/store/user/jhgoh/CATTools/sync/v7-6-3/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',]
#process.source.fileNames = ['file:../../../catdata_20160315/catTuple.root']
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TTTo2L2Nu_13TeV-powheg/v8-0-0_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/160705_215520/0000/catTuple_1.root']


from CATTools.CatAnalyzer.leptonSF_cff import *

process.cattree = cms.EDAnalyzer("JetChargeAnalyzer",
    muon = cms.PSet(
        src = cms.InputTag("catMuons"),
    ),
    electron = cms.PSet(
        src = cms.InputTag("catElectrons"),
    ),
    jets = cms.InputTag("catJets"),
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("jet_charge_tree.root"
))

process.p = cms.Path(process.cattree)
if ( process.maxEvents.input <0 or process.maxEvents > 5000) :
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = True
