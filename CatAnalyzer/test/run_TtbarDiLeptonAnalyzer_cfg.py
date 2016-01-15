import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151127_200613/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/MuonEG/v7-4-6_Run2015C_25ns-05Oct2015-v1/151127_195051/0000/catTuple_1.root')

useSilver = False
catmet = 'catMETs'
lumiMask = 'lumiMask'
pileupWeight = 'pileupWeight'
if useSilver:
    catmet = 'catMETsNoHF'
    lumiMask = 'lumiMaskSilver'
    pileupWeight = 'pileupWeightSilver'
    
process.filterRECO = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::PAT"),
    secondaryTriggerResults = cms.InputTag("TriggerResults::RECO"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("and"),
    triggersToMatch = cms.vstring(
        "CSCTightHaloFilter",
        #"EcalDeadCellTriggerPrimitiveFilter",
        #"HBHENoiseFilter",
        "eeBadScFilter",
        "goodVertices",
    ),
    doFilter = cms.bool(False),
)

process.filterTrigMUEL = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("or"),
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
    ),
    doFilter = cms.bool(False),
)

process.filterTrigELEL = process.filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
 ## "HLT_Ele23_WPLoose_Gsf_v",
    ),
)

process.filterTrigMUMU = process.filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
 ## "HLT_IsoMu20_v",
 ## "HLT_IsoTkMu20_v",
 ## "HLT_IsoMu27_v",
    ),
)

process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff")
from CATTools.CatAnalyzer.leptonSF_cff import *

process.cattree = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag(lumiMask),
    genweight = cms.InputTag("genWeight","genWeight"),
    pdfweight = cms.InputTag("genWeight","pdfWeights"),
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
        effSF = electronSFWP90,
    ),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag(catmet),
    mcLabel = cms.InputTag("prunedGenParticles"),
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    pseudoTop = cms.InputTag("pseudoTop"),
    
    solver = process.ttbarDileptonKinAlgoPSetCMSKin,
    #solver = process.ttbarDileptonKinAlgoPSetDESYSmeared,
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
