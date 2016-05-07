import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-3/MuonEG_Run2015D-16Dec2015-v1.root',]
#process.source.fileNames = ['file:/xrootd/store/user/jhgoh/CATTools/sync/v7-6-3/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',]
process.source.fileNames = ['file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_1.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_2.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_3.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_4.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_5.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_6.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_7.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_8.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_9.root',
'file:/xrootd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-6-3_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/160221_150606/0000/catTuple_10.root',
]

useSilver = False
catmet = 'catMETs'
lumiMask = 'lumiMask'
pileupWeight = 'pileupWeight'
if useSilver:
    catmet = 'catMETsNoHF'
    lumiMask = 'lumiMaskSilver'
    pileupWeight = 'pileupWeightSilver'

process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff")
process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.topPtWeightProducer_cfi")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
from CATTools.CatAnalyzer.leptonSF_cff import *

process.ttbarDileptonKinAlgoPSetDESYSmeared.inputTemplatePath = cms.string("CATTools/CatAnalyzer/data/desyKinRecoInput.root")
process.ttbarDileptonKinAlgoPSetDESYSmeared.maxLBMass = cms.double(180)
process.ttbarDileptonKinAlgoPSetDESYSmeared.mTopInput = cms.double(172.5)
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop = process.ttbarDileptonKinAlgoPSetDESYSmeared.clone()
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop.inputTemplatePath = cms.string("CATTools/CatAnalyzer/data/KoreaKinRecoInput_pseudo.root")
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop.maxLBMass = cms.double(360)
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop.mTopInput = cms.double(172.5)

process.cattree = cms.EDAnalyzer("dileptonCommon",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag(lumiMask),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweights = cms.InputTag("flatGenWeights","pdf"),		
    scaleupweights = cms.InputTag("flatGenWeights","scaleup"),
    scaledownweights = cms.InputTag("flatGenWeights","scaledown"),
    topPtWeight = cms.InputTag("topPtWeight"),
    puweight = cms.InputTag(pileupWeight),
    puweight_up = cms.InputTag(pileupWeight,"up"),
    puweight_dn = cms.InputTag(pileupWeight,"dn"),
    trigMUEL = cms.InputTag("filterTrigMUEL"),
    trigMUMU = cms.InputTag("filterTrigMUMU"),
    trigELEL = cms.InputTag("filterTrigELEL"),

    vertices = cms.InputTag("catVertex"),
    muon = cms.PSet(
        src = cms.InputTag("catMuons"),
        effSF = muonSFTight,
    ),
    electron = cms.PSet(
        src = cms.InputTag("catElectrons"),
        effSF = electronSFCutBasedIDMediumWP,#electronSFWP90,
    ),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag(catmet),
    mcLabel = cms.InputTag("prunedGenParticles"),
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    pseudoTop = cms.InputTag("pseudoTop"),
    
    #solver = process.ttbarDileptonKinAlgoPSetCMSKin,
    solver = process.ttbarDileptonKinAlgoPSetDESYSmeared,
    solverPseudoTop = process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop,
    #solver = process.ttbarDileptonKinAlgoPSetDESYMassLoop,
)
#process.cattree.solver.tMassStep = 1
if cms.string('DESYSmeared') == process.cattree.solver.algo:
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        cattree = cms.PSet(
            initialSeed = cms.untracked.uint32(123456),
            engineName = cms.untracked.string('TRandom3')
        )
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root"
))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
