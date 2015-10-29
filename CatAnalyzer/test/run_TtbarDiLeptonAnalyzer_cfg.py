import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-4_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_143547/0000/catTuple_96.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/DoubleMuon/v7-4-4_Run2015C_25ns-05Oct2015-v1/151023_165157/0000/catTuple_10.root')

process.filterRECO = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::PAT"),
    secondaryTriggerResults = cms.InputTag("TriggerResults::RECO"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("and"),
    triggersToMatch = cms.vstring(
        "CSCTightHaloFilter",
        "EcalDeadCellTriggerPrimitiveFilter",
        "HBHENoiseFilter",
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
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    ),
)

process.filterTrigMUMU = process.filterTrigMUEL.clone(
    triggersToMatch = cms.vstring(
      "HLT_Mu17_Mu8_DZ_v",
      "HLT_Mu17_TkMu8_DZ_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
    ),
)

process.load("CATTools.CatAnalyzer.ttbarDileptonKinSolutionAlgos_cff")

process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    puweight = cms.InputTag("pileupWeight"),
    trigMUEL = cms.InputTag("filterTrigMUEL"),
    trigMUMU = cms.InputTag("filterTrigMUMU"),
    trigELEL = cms.InputTag("filterTrigELEL"),

    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    #mets = cms.InputTag("catMETsNoHF"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    pseudoTop = cms.InputTag("pseudoTop"),
    
    solver = process.ttbarDileptonKinAlgoPSetCMSKin,
)
process.ttll.solver.tMassStep = 1
#process.ttll.solver.algo = cms.string('DESYMASSLOOP')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("top.root"
))

process.p = cms.Path(process.ttll)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
