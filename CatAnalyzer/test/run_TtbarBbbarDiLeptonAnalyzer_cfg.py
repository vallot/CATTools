import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-1_RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160123_212836/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/151222_133640/0000/catTuple_1.root')

catmet = 'catMETs' #NoHF'
lumiMask = 'lumiMask'

process.load("CATTools.CatAnalyzer.filters_cff")
from CATTools.CatAnalyzer.leptonSF_cff import *

process.cattree = cms.EDAnalyzer("TtbarBbbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    genweight = cms.InputTag("genWeight","genWeight"),
    genweightQ = cms.InputTag("genWeight","Q"),
    genweightPDF = cms.InputTag("genWeight","pdfWeights"),

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

    elecSF = electronSFWP90,
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
    #solver = process.ttbarDileptonKinAlgoPSetCMSKin,
    #solver = process.ttbarDileptonKinAlgoPSetDESYSmeared,
    #solver = process.ttbarDileptonKinAlgoPSetDESYMassLoop,
)
#process.cattree.solver.tMassStep = 1
#if cms.string('DESYSmeared') == process.cattree.solver.algo:
#    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#        cattree = cms.PSet(
#            initialSeed = cms.untracked.uint32(123456),
#            engineName = cms.untracked.string('TRandom3')
#        )
#    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root"
))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000


