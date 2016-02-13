#------------------------------------------------------------------
#------------------------------------------------------------------
# Data or MC Sample
runOnMC      = True
# runOnTTbarMC == 0, No ttbar
# runOnTTbarMC == 1, ttbar Signal
# runOnTTbarMC == 2, ttbar Background
runOnTTbarMC = 1
#------------------------------------------------------------------
#------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
process = cms.Process("ttbbLepJets")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('ttbbLepJets')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:ttbar_PowhegPythia.root',
        #'file:catTuple_1.root',
        #'/store/user/jhgoh/CATTools/sync/v7-6-2/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',
    )
)

# PUReWeight
# process.load("CATTools.CatProducer.pileupWeight_cff")
# from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
# process.pileupWeight.weightingMethod = "RedoWeight"
# process.pileupWeight.pileupRD = pileupWeightMap["Run2015_25nsV1"]
# process.pileupWeight.pileupUp = pileupWeightMap["Run2015Up_25nsV1"]
# process.pileupWeight.pileupDn = pileupWeightMap["Run2015Dn_25nsV1"]

# json file (Only Data)
# ReReco JSON file taken from: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt
# import FWCore.PythonUtilities.LumiList as LumiList
# process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

# Lepton SF
from CATTools.CatAnalyzer.leptonSF_cff import *

process.ttbbLepJets = cms.EDAnalyzer('ttbbLepJetsAnalyzer',
                                     sampleLabel       = cms.untracked.bool(runOnMC),
                                     TTbarSampleLabel  = cms.untracked.int32(runOnTTbarMC),
                                     genWeightLabel    = cms.InputTag("genWeight:genWeight"),
                                     ScaleWeightLabel  = cms.InputTag("genWeight:scaleWeights"),
                                     genLabel          = cms.InputTag("prunedGenParticles"),
                                     genJetLabel       = cms.InputTag("slimmedGenJets"),
                                     genttbarCatLabel  = cms.InputTag("GenTtbarCategories:genTtbarId"),
                                     genttbarConeCatLabel  = cms.InputTag("catGenTops"),
                                     muonLabel         = cms.InputTag("catMuons"),
                                     muonSF = muonSFTight,
                                     electronLabel     = cms.InputTag("catElectrons"),
                                     elecSF = electronSFWP90,
                                     jetLabel          = cms.InputTag("catJets"),
                                     metLabel          = cms.InputTag("catMETs"),
                                     pvLabel           = cms.InputTag("catVertex:nGoodPV"),
                                     puWeightLabel     = cms.InputTag("pileupWeight"),
                                     puUpWeightLabel   = cms.InputTag("pileupWeight:up"),
                                     puDownWeightLabel = cms.InputTag("pileupWeight:dn"),
                                     triggerBits       = cms.InputTag("TriggerResults::HLT"), 
                                     triggerObjects    = cms.InputTag("catTrigger"), 
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Tree_ttbbLepJets.root')
                                   )


# process.Tracer = cms.Service("Tracer")
# process.dump=cms.EDAnalyzer('EventContentAnalyzer')
# process.p = cms.Path(process.demo*process.dump)
# process.p = cms.Path(process.pileupWeight*
#                      process.ttbarSingleLepton)
process.p = cms.Path(process.ttbbLepJets)
