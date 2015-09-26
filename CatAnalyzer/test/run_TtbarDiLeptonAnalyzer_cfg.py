import FWCore.ParameterSet.Config as cms

process = cms.Process("TtbarDiLeptonAnalyzer")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

datadir = '/xrootd/store/group/CAT/MuonEG/v7-3-6_Run2015B-PromptReco-v1/150922_133849/0000/'
datadir = '/xrootd/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/'
datadir = '/xrootd/store/group/CAT/DoubleMuon/v7-3-6_Run2015B-PromptReco-v1/150922_133736/0000/'

import os
for f in os.listdir(datadir):
    if ".root" in f:
        process.source.fileNames.append("file:"+datadir+f)

#process.source.fileNames.append('file:/cms/scratch/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_1.root')
#process.source.fileNames.append('file:/cms/scratch/CAT/WW_TuneCUETP8M1_13TeV-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150805_203816/0000/catTuple_1.root')
#process.source.fileNames.append('file:/afs/cern.ch/user/j/jlee/cat74/src/CATTools/CatProducer/prod/catTuple.root')
#process.source.fileNames.append('/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-4_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150810_215031/0000/catTuple_101.root')
#process.source.fileNames.append('file:/afs/cern.ch/user/j/jlee/test/cat74/src/CATTools/CatProducer/prod/catTuple.root')

#lumiFile = 'Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt'
lumiFile = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'

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

process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    goodVertices = cms.InputTag("catTrigger", "goodVertices"),
    CSCTightHaloFilter = cms.InputTag("catTrigger", "CSCTightHaloFilter"),
    HBHENoiseFilter = cms.InputTag("catTrigger", "HBHENoiseFilter"),
    eeBadScFilter = cms.InputTag("catTrigger", "eeBadScFilter"),
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerObjects = cms.InputTag("catTrigger"),

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

    pseudoTop_jets = cms.InputTag("pseudoTop","jets"),
    pseudoTop_leptons = cms.InputTag("pseudoTop","leptons"),
    pseudoTop = cms.InputTag("pseudoTop"),
    pseudoTop_neutrinos = cms.InputTag("pseudoTop","neutrinos"),
    pseudoTop_mets = cms.InputTag("pseudoTop","mets"),
    
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
