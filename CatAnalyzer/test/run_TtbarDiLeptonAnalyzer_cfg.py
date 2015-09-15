import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("TtbarDiLeptonAnalyzer")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

"""
###### discard this line when use for data sample ##########
process.partonTop = cms.EDProducer("PartonTopProducer",
    genParticles = cms.InputTag("prunedGenParticles"),
    jetMinPt = cms.double(20),
    jetMaxEta = cms.double(2.5),
    jetConeSize = cms.double(0.4),
)
###### discard this line when use for data sample ###########
"""

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

#process.source.fileNames.append('file:/cms/scratch/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_1.root')
#process.source.fileNames.append('file:/cms/scratch/CAT/WW_TuneCUETP8M1_13TeV-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150805_203816/0000/catTuple_1.root')
#process.source.fileNames.append('file:/afs/cern.ch/user/j/jlee/cat74/src/CATTools/CatProducer/prod/catTuple.root')
#process.source.fileNames.append('/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150820_215807/0000/catTuple_193.root')
process.source.fileNames.append('/store/group/CAT/MuonEG/desySync_Run2015B-17Jul2015-v1/150905_193954/0000/catTuple_1.root')

#process.source.lumisToProcess = LumiList.LumiList(filename = 'rereco_JSON.txt').getVLuminosityBlockRange()
process.source.lumisToProcess = LumiList.LumiList(filename = 'prompt_JSON.txt').getVLuminosityBlockRange()


process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    goodVertices = cms.InputTag("catTrigger", "goodVertices"),
    CSCTightHaloFilter = cms.InputTag("catTrigger", "CSCTightHaloFilter"),
    HBHENoiseFilter = cms.InputTag("catTrigger", "HBHENoiseFilter"),
    eeBadScFilter = cms.InputTag("catTrigger", "eeBadScFilter"),
    HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL = cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL"),
    HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL = cms.InputTag("catTrigger", "HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL"),

    vertices = cms.InputTag("catVertex"),
    #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    noHFmets = cms.InputTag("catMETsNoHF"),
    mcLabel = cms.InputTag("prunedGenParticles"),
    triggers = cms.InputTag("catTrigger"),
    
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
