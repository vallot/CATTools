import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-1_RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160203_113158/0000/catTuple_1.root')
process.source.fileNames.append('/store/user/jhgoh/CATTools/sync/v7-6-1/TTbarXSecSynchronization_76X_MC_TT_powheg.root')
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-1/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-1/DoubleEG_Run2015D-16Dec2015-v2.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-1/DoubleMuon_Run2015D-16Dec2015-v1.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-1/MuonEG_Run2015D-16Dec2015-v1.root',]

#import os
#useGold = True
#isRunData = False
catmet = 'catMETs' #NoHF'
lumiMask = 'lumiMask'
#if useGold:
#    catmet = 'catMETs'
#    if isRunData:
#        #lumiFile = 'Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#        lumiFile = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#        from FWCore.PythonUtilities.LumiList import LumiList
#        lumiList = LumiList(os.environ["CMSSW_BASE"]+'/src/CATTools/CatProducer/prod/LumiMask/'+lumiFile)
#        process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()
    
#    process.load("CATTools.CatProducer.pileupWeight_cff")
#    from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
#    process.pileupWeight.weightingMethod = "RedoWeight"
#    process.pileupWeight.pileupRD = pileupWeightMap["Run2015_25nsV1"]
#    process.pileupWeight.pileupUp = pileupWeightMap["Run2015Up_25nsV1"]
#    process.pileupWeight.pileupDn = pileupWeightMap["Run2015Dn_25nsV1"]

#process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff")

####for running genTop on the fly. however it is running slowly.
#process.load("CATTools.CatProducer.genTopProducer_cfi")
#from CATTools.CatProducer.Tools.tools import genHFTool
#genHFTool(process,True)

process.load("CATTools.CatAnalyzer.filters_cff")

##for only ttbar signal mc sample
process.load("CATTools.CatAnalyzer.topPtWeightProducer_cfi")

from CATTools.CatAnalyzer.leptonSF_cff import *

process.cattree = cms.EDAnalyzer("TtbarBbbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    genweight = cms.InputTag("genWeight","genWeight"),
    genweightQ = cms.InputTag("genWeight","Q"),
    #genweightPDF = cms.InputTag("genWeight","pdfWeights"),
    genweightPDF = cms.InputTag("genWeight","otherWeights"),
    scaleweight = cms.InputTag("genWeight","scaleWeights"),
    topPtWeight = cms.InputTag("topPtWeight"),

    lumiSelection = cms.InputTag(lumiMask),
    puweight = cms.InputTag("pileupWeight"),
    puweightUp = cms.InputTag("pileupWeight","up"),
    puweightDown = cms.InputTag("pileupWeight","dn"),
    trigMUEL = cms.InputTag("filterTrigMUEL"),
    trigMUMU = cms.InputTag("filterTrigMUMU"),
    trigELEL = cms.InputTag("filterTrigELEL"),

    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag(catmet),
    mcLabel = cms.InputTag("prunedGenParticles"),

    #elecSF = electronSFWP90,
    elecSF = electronSFCutBasedIDMediumWP,
    muonSF = muonSFTight,
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    pseudoTop = cms.InputTag("pseudoTop"),

    genTtbarId = cms.InputTag("GenTtbarCategories", "genTtbarId"),
    genTtbarId30 = cms.InputTag("GenTtbarCategories30", "genTtbarId"),
    genTtbarId40 = cms.InputTag("GenTtbarCategories40", "genTtbarId"),

    GenJets = cms.InputTag("slimmedGenJets"),
    GenParticles = cms.InputTag("prunedGenParticles"),
    GenTop = cms.InputTag("catGenTops"),
    #GenTop = cms.InputTag("catGenTops","","TtbarDiLeptonAnalyzer"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root"
))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000


