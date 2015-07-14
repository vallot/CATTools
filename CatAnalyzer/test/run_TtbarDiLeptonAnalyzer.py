import FWCore.ParameterSet.Config as cms

process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#for i in xrange(1,101):
#    process.source.fileNames.append('file:/cms/data/xrd/store/user/jlee/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/cat74v2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150713_164609/0000/catTuple_%d.root' % i)

process.source.fileNames.append('file:/cms/data/xrd/store/user/jlee/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/cat74v2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150713_164609/0000/catTuple_1.root')

process.ttll = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("top.root"
))

process.p = cms.Path(process.ttll)
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
