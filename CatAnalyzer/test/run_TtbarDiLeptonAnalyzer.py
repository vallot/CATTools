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

process.source.fileNames.append('/store/group/CAT/SingleMuon/v7-3-6_Run2015B-17Jul2015-v1/150820_215426/0000/catTuple_6.root')
#process.source.fileNames.append('/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150820_215807/0000/catTuple_193.root')

process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt').getVLuminosityBlockRange()

process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    vertices = cms.InputTag("catVertex"),
    #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
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
