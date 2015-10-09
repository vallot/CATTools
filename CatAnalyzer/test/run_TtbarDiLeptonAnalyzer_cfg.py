import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("TtbarDiLeptonAnalyzer")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

#datadir = '/xrootd/store/group/CAT/MuonEG/v7-3-6_Run2015B-PromptReco-v1/150922_133849/0000/'
#datadir = '/xrootd/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/'
#datadir = '/xrootd/store/group/CAT/DoubleMuon/v7-3-6_Run2015B-PromptReco-v1/150922_133736/0000/'
#
#for f in os.listdir(datadir):
#    if ".root" in f:
#        process.source.fileNames.append("file:"+datadir+f)
#

#process.source.fileNames = ['file:catTuple.root']
process.source.fileNames.append('/store/group/CAT/MuonEG/v7-4-2_Run2015C-PromptReco-v1/150923_202331/0000/catTuple_2.root')
process.source.fileNames.append('/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-2_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150923_215647/0001/catTuple_1046.root')
process.source.fileNames.append('/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-2_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150923_215647/0001/catTuple_1047.root')
process.source.fileNames.append('/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-2_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150923_215647/0001/catTuple_1048.root')
process.source.fileNames.append('/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-2_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150923_215647/0001/catTuple_1049.root')

#lumiFile = 'Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt'
lumiFile = 'Cert_246908-257599_13TeV_PromptReco_Collisions15_25ns_JSON.txt'

runOnMC = True
for i in process.source.fileNames:
    if 'Run2015' in i:
        runOnMC=False
if not runOnMC:
    from FWCore.PythonUtilities.LumiList import LumiList
    lumiList = LumiList(os.environ["CMSSW_BASE"]+'/src/CATTools/CatProducer/prod/LumiMask/'+lumiFile)    
    process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()    
    
if runOnMC:
    process.partonTop = cms.EDProducer("PartonTopProducer",
        genParticles = cms.InputTag("prunedGenParticles"),
        jetMinPt = cms.double(20),
        jetMaxEta = cms.double(2.5),
        jetConeSize = cms.double(0.4),
    )

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

process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    trigMUEL = cms.InputTag("filterTrigMUEL"),
    trigMUMU = cms.InputTag("filterTrigMUMU"),
    trigELEL = cms.InputTag("filterTrigELEL"),

    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    #mets = cms.InputTag("catMETs"),
    mets = cms.InputTag("catMETsNoHF"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    
    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    #isTTbarMC = cms.bool(True),
    isTTbarMC = cms.bool(False),
    pseudoTop = cms.InputTag("pseudoTop"),
    
    tmassbegin = cms.double(100),
    tmassend   = cms.double(300),
    tmassstep  = cms.double(  1),
    neutrino_parameters = cms.vdouble(27.23,53.88,19.92,53.89,19.9)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("top.root"
))

process.p = cms.Path(process.ttll)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
