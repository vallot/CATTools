import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v8-0-2_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/161104_204003/0002/catTuple_2111.root'] # TT_powheg v8-0-2
#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v8-0-1_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/160822_144100/0002/catTuple_2112.root'] # TT_powheg v8-0-1
#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/ZZ_TuneCUETP8M1_13TeV-pythia8/v8-0-1_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160810_124436/0000/catTuple_22.root'] # ZZ v8-0-1
#process.source.fileNames = ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-1_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160810_120220/0000/catTuple_85.root'] # DYJets v8-0-1

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

process.agen = cms.EDAnalyzer("CATGenTopAnalysis",
    weightIndex = cms.int32(-1),
    weight = cms.InputTag("flatGenWeights"),
    channel = cms.InputTag("partonTop","channel"),
    modes = cms.InputTag("partonTop", "modes"),
    partonTop = cms.InputTag("partonTop"),
    pseudoTop = cms.InputTag("pseudoTop"),
    filterTaus = cms.bool(False),
)

process.cattree = cms.EDAnalyzer("CATDstarAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag("lumiMask"),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweight = cms.InputTag("flatGenWeights","pdf"),
    scaleupweight = cms.InputTag("flatGenWeights","scaleup"),
    scaledownweight = cms.InputTag("flatGenWeights","scaledown"),
    topPtWeight = cms.InputTag("topPtWeight"),
    puweight = cms.InputTag("pileupWeight"),
    puweight_up = cms.InputTag("pileupWeight","up"),
    puweight_dn = cms.InputTag("pileupWeight","dn"),
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
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    # input collection
    d0s    = cms.InputTag("catDstars","D0Cand"),
    dstars = cms.InputTag("catDstars","DstarCand"),
    Jpsis  = cms.InputTag("catDstars","JpsiCand"),
    matchingDeltaR = cms.double(0.15),
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
    fileName = cms.string("cattree_CATDstar.root"
))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
process.options.wantSummary = True
