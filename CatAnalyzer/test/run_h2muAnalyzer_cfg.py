import FWCore.ParameterSet.Config as cms
process = cms.Process("h2muAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
from CATTools.Validation.commonTestInput_cff import commonTestCATTuples
#process.source.fileNames = commonTestCATTuples["sig"]
process.source.fileNames = cms.untracked.vstring(
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_1.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_2.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_3.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_4.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_5.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_6.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_7.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_8.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_9.root",
    "file:/xrootd/store/group/CAT/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_114412/0000/catTuple_10.root" )

#txtfile = '../data/dataset/dataset_SingleMuon_Run2016E.txt'
#txtfile = '../data/dataset/dataset_DYJets.txt'
#f = open(txtfile)
#for line in f:
#    if '#' not in line:
#        process.source.fileNames.append(line)
#print process.source.fileNames

process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
from CATTools.CatAnalyzer.leptonSF_cff import *

process.cattree = cms.EDAnalyzer("h2muAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag("lumiMask"),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweight = cms.InputTag("flatGenWeights","pdf"),
    scaleweight = cms.InputTag("flatGenWeights","scaleWeights"),
    puweight = cms.InputTag('pileupWeight'),
    puweight_up = cms.InputTag('pileupWeight',"up"),
    puweight_dn = cms.InputTag('pileupWeight',"dn"),
    vertices = cms.InputTag("catVertex"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag('catMETs'),
    mcLabel = cms.InputTag("prunedGenParticles"),
    triggerBits = cms.VInputTag(
        cms.InputTag("TriggerResults","","HLT2"),# due to reHLT, this is the first choice 
        cms.InputTag("TriggerResults","","HLT"),# if above is not found, falls to default 
    ),
    triggerObjects = cms.InputTag("catTrigger"),
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
    fileName = cms.string("cattree.root"))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
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
