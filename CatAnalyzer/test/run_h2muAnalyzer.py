import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('dataset', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "dataset: 1  default")
options.parseArguments()

process = cms.Process("h2muAnalyzer")
savename ="h2mu"
datadir ="/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-0_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150720_065744/0000/"

if options.dataset == 1:
    savename ="DYJetsToLL_M-50"
    datadir ="/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-0_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150720_065744/0000/"
if options.dataset == 2:
    savename ="TTJets_TuneCUETP8M1_13TeV-madgraphMLM"
    datadir ="/store/group/CAT/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/v7-3-0_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150720_065809/0000/"
if options.dataset == 3:
    savename ="ZZ_TuneCUETP8M1_13TeV"
    datadir ="/store/group/CAT/ZZ_TuneCUETP8M1_13TeV-pythia8/v7-3-0_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150720_070032/0000/"
if options.dataset == 4:
    savename ="WZ_TuneCUETP8M1_13TeV"
    datadir ="/store/group/CAT/WZ_TuneCUETP8M1_13TeV-pythia8/v7-3-0_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150720_070007/0000/"
if options.dataset == 5:
    savename ="WW_TuneCUETP8M1_13TeV"
    datadir ="/store/group/CAT/WW_TuneCUETP8M1_13TeV-pythia8/v7-3-0_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150720_065935/0000/"
        
if options.dataset == 6:
    savename ="SingleMuon"
    datadir ="/store/group/CAT/SingleMuon/v7-3-0_Run2015B-PromptReco-v1/150720_060727/0000/"
if options.dataset == 7:
    savename ="DoubleMuon"
    datadir ="/store/group/CAT/DoubleMuon/v7-3-0_Run2015B-PromptReco-v1/150720_060849/0000"

savename+=".root"

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
import os
for f in os.listdir(datadir):
    process.source.fileNames.append("file:"+datadir+f)

process.h2mu = cms.EDAnalyzer("h2muAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(savename
))

process.p = cms.Path(process.h2mu)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
