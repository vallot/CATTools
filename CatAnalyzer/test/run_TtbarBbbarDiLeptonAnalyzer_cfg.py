import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151127_200613/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd//store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151130_231748/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd//store/group/CAT/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151127_200111/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151130_231748/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/DoubleMuon/v7-4-6_Run2015D-PromptReco-v4/151127_194953/0000/catTuple_1.root')
#process.source.fileNames.append('root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/DoubleMuon/v7-4-4_Run2015C_25ns-05Oct2015-v1/151023_165157/0000/catTuple_10.root')

#import os
#useGold = True
#isRunData = False
catmet = 'catMETsNoHF'
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
process.filterRECO = cms.EDFilter("CATTriggerBitCombiner",
    triggerResults = cms.InputTag("TriggerResults::PAT"),
    secondaryTriggerResults = cms.InputTag("TriggerResults::RECO"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineBy = cms.string("and"),
    triggersToMatch = cms.vstring(
        "CSCTightHaloFilter",
#        "EcalDeadCellTriggerPrimitiveFilter",
#        "HBHENoiseFilter",
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
#      "HLT_Mu17_Mu8_DZ_v",
#      "HLT_Mu17_TkMu8_DZ_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
    ),
)

#process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff")

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
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    pseudoTop = cms.InputTag("pseudoTop"),

    genTtbarId = cms.InputTag("GenTtbarCategories", "genTtbarId"),
    genTtbarId30 = cms.InputTag("GenTtbarCategories30", "genTtbarId"),
    genTtbarId40 = cms.InputTag("GenTtbarCategories40", "genTtbarId"),

    GenJets = cms.InputTag("slimmedGenJets"),
    GenParticles = cms.InputTag("prunedGenParticles"),
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


