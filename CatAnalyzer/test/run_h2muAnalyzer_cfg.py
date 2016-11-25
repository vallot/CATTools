import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("h2muAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),)
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/DoubleEG_Run2015D-16Dec2015-v2.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/DoubleMuon_Run2015D-16Dec2015-v1.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/MuonEG_Run2015D-16Dec2015-v1.root',]
#process.source.fileNames = ['file:%s/src/CATTools/CatProducer/prod/catTuple.root'%os.environ["CMSSW_BASE"]]

process.source.fileNames = []
for i in range(1):
  process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/SingleMuon/v8-0-1_Run2016D-PromptReco-v2/160810_112606/0000/catTuple_%s.root'%(i+1))
#process.source.fileNames=['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/SingleMuon/v8-0-1_Run2016D-PromptReco-v2/160810_112606/0000/catTuple_1.root']


#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/VBF_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-1_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v2/160810_125753/0000/catTuple_1.root']
#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-1_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160810_120220/0000/catTuple_324.root',]
#txtfile = '../data/dataset/dataset_SingleMuon_Run2016E.txt'
#txtfile = '../data/dataset/dataset_DYJets.txt'
#f = open(txtfile)
#for line in f:
#    if '#' not in line:
#        process.source.fileNames.append(line)
#print process.source.fileNames
    
catmet = 'catMETs'
lumiMask = 'lumiMask'
pileupWeight = 'pileupWeight'

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
from CATTools.CatAnalyzer.leptonSF_cff import *

process.cattree = cms.EDAnalyzer("h2muAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag(lumiMask),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweight = cms.InputTag("flatGenWeights","pdfWeights"),
    scaleweight = cms.InputTag("flatGenWeights","scaleWeights"),
    puweight = cms.InputTag(pileupWeight),
    puweight_up = cms.InputTag(pileupWeight,"up"),
    puweight_dn = cms.InputTag(pileupWeight,"dn"),
    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag(catmet),
    mcLabel = cms.InputTag("prunedGenParticles"),
    #triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerBits = cms.VInputTag(
        cms.InputTag("TriggerResults","","HLT2"),# due to reHLT, this is the first choice 
        cms.InputTag("TriggerResults","","HLT"),# if above is not found, falls to default 
    ),
    triggerObjects = cms.InputTag("catTrigger"),
    #triggerObjects = cms.InputTag("selectedPatTrigger"),
    muon = cms.PSet(
        src = cms.InputTag("catMuons"),
        effSF = muonSFTight,
    ),
    electron = cms.PSet(
        src = cms.InputTag("catElectrons"),
        effSF = electronSFCutBasedIDMediumWP,#electronSFWP90,
    ),    
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root")
)

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
"""
process.MessageLogger = cms.Service("MessageLogger",
    destinations   = cms.untracked.vstring(
        'detailedInfo' 
    ),
    detailedInfo   = cms.untracked.PSet(
        #reportEvery = cms.untracked.int32(50000),
        extension = cms.untracked.string('.txt') 
    )
)
"""
