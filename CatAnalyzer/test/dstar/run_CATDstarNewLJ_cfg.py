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
process.source.fileNames = ['file:catTuple_LJ.root']


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

## Redo the pileup weight - necessary for v765 production
process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeight = process.pileupWeight.clone()
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.redoPileupWeight.weightingMethod = "RedoWeight"
process.redoPileupWeight.pileupMC = pileupWeightMap["2015_25ns_FallMC"]
process.redoPileupWeight.pileupRD = pileupWeightMap["Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON"]
process.redoPileupWeight.pileupUp = pileupWeightMap["Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Up"]
process.redoPileupWeight.pileupDn = pileupWeightMap["Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Dn"]
pileupWeight = 'redoPileupWeight'

process.cattree = cms.EDAnalyzer("CATDstarSemiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag(lumiMask),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweight = cms.InputTag("flatGenWeights","pdf"),		
    scaleupweight = cms.InputTag("flatGenWeights","scaleup"),
    scaledownweight = cms.InputTag("flatGenWeights","scaledown"),
    topPtWeight = cms.InputTag("topPtWeight"),
    puweight = cms.InputTag(pileupWeight),
    puweight_up = cms.InputTag(pileupWeight,"up"),
    puweight_dn = cms.InputTag(pileupWeight,"dn"),
    trigMUJET = cms.InputTag("filterTrigMUJET"),
    trigELJET = cms.InputTag("filterTrigELJET"),

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
    # input collection
    d0s    = cms.InputTag("catDstars","D0Cand"),
    dstars = cms.InputTag("catDstars","DstarCand"),
    matchingDeltaR = cms.double(0.5),
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    pseudoTop = cms.InputTag("pseudoTop"),

    MuonPtCut = cms.double(26.),
    MuonEtaCut = cms.double(2.1),
    MuonIDCut = cms.string("tight"), # tight or loose
    MuonIsoCut = cms.double(0.15),
    vetoMuonPtCut = cms.double(15.),
    vetoMuonEtaCut = cms.double(2.4),
    vetoMuonIDCut = cms.string("loose"),
    vetoMuonIsoCut = cms.double(0.25),
    ElectronPtCut = cms.double(30.),
    ElectronEtaCut = cms.double(2.4),
    ElectronIDCut = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    ElectronIsoCut = cms.double(999),
    vetoElectronPtCut = cms.double(15.),
    vetoElectronEtaCut = cms.double(2.4),
    vetoElectronIDCut = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    vetoElectronIsoCut = cms.double(999.),
    JetPtCut = cms.double(30.),
    JetEtaCut = cms.double(2.4),
    JetCheckOverlap = cms.bool(True),

    bJetCSV = cms.string("medium"),
    

    
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree_CATDstarNewLJ.root"
))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
process.options.wantSummary = True
