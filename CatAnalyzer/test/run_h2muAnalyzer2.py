import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('dataset', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "dataset: 1  default")
options.parseArguments()

process = cms.Process("h2muAnalyzer")
savename ="h2mu"
datadir = "/xrootd/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150820_215635/0000/"
first_dir = '/xrootd/store/group/CAT/'

if options.dataset == 1:
    savename ="DYJetsToLL_M-50"
    datadir = first_dir+'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150820_215635/0000/catTuple_100.root'
if options.dataset == 2:
    savename ="TTJets_TuneCUETP8M1_13TeV-madgraphMLM"
    datadir = first_dir+'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150820_215713/0000/catTuple_100.root'
if options.dataset == 3:
    savename ="ZZ_TuneCUETP8M1_13TeV"
    datadir = first_dir+'ZZ_TuneCUETP8M1_13TeV-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150820_215925/0000/catTuple_20.root'
if options.dataset == 4:
    savename ="WZ_TuneCUETP8M1_13TeV"
    datadir = first_dir+'WZ_TuneCUETP8M1_13TeV-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150820_215857/0000/catTuple_35.root'
if options.dataset == 5:
    savename ="WW_TuneCUETP8M1_13TeV"
    datadir = first_dir+'WW_TuneCUETP8M1_13TeV-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150820_215830/0000/catTuple_39.root'
if options.dataset == 6:
    savename ="SingleMuon"
    datadir = first_dir+'SingleMuon/v7-3-6_Run2015B-PromptReco-v1/150820_215216/0000/catTuple_31.root'
if options.dataset == 7:
    savename ="DoubleMuon"
    datadir = first_dir+'DoubleMuon/v7-3-6_Run2015B-PromptReco-v1/150820_215326/0000/catTuple_22.root'

savename+=".root"

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

#import os
#for f in os.listdir(datadir):
process.source.fileNames.append("file:"+datadir)

process.h2mu = cms.EDAnalyzer("h2muAnalyzer",
    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(savename)
)

process.p = cms.Path(process.h2mu)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
