import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('dataset', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "dataset: 1  default")
options.parseArguments()

process = cms.Process("h2muAnalyzer")
savename ="h2mu"
datadir ="/cms/scratch/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150805_203735/0000/"
if options.dataset == 1:
    savename ="DYJetsToLL_M-50"
    datadir ="/cms/scratch/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150805_203735/0000/"
if options.dataset == 2:
    savename ="TTJets_TuneCUETP8M1_13TeV-madgraphMLM"
    datadir ="/cms/scratch/CAT/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150805_203745/0000/"
if options.dataset == 3:
    savename ="ZZ_TuneCUETP8M1_13TeV"
    datadir ="/cms/scratch/CAT/ZZ_TuneCUETP8M1_13TeV-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150805_203839/0000/"
if options.dataset == 4:
    savename ="WZ_TuneCUETP8M1_13TeV"
    datadir ="/cms/scratch/CAT/WZ_TuneCUETP8M1_13TeV-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/150805_203826/0000/"
if options.dataset == 5:
    savename ="WW_TuneCUETP8M1_13TeV"
    datadir ="/cms/scratch/CAT/WW_TuneCUETP8M1_13TeV-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150805_203816/0000/"
        
if options.dataset == 6:
    savename ="SingleMuon"
    datadir ="/cms/scratch/CAT/SingleMuon/v7-3-2_Run2015B-PromptReco-v1/150806_164946/0000/"
if options.dataset == 7:
    savename ="DoubleMuon"
    datadir ="/cms/scratch/CAT/DoubleMuon/v7-3-2_Run2015B-PromptReco-v1/150806_165018/0000/"

savename+=".root"

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
import os
for f in os.listdir(datadir):
    process.source.fileNames.append("file:"+datadir+f)

process.h2mu = cms.EDAnalyzer("h2muAnalyzer",
    vertices = cms.InputTag("catVertex"),
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
