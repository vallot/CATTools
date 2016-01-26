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
        #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151127_200111/0000/catTuple_102.root'
        'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-6_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151222_133941/0000/catTuple_1.root'
        #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/SingleMuon/v7-4-5_Run2015D-05Oct2015-v1/151109_230333/0000/catTuple_1.root'
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
# import FWCore.PythonUtilities.LumiList as LumiList
# process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

process.ttbbLepJets = cms.EDAnalyzer('ttbbLepJetsAnalyzer',
                                     sampleLabel       = cms.untracked.bool(runOnMC),
                                     TTbarSampleLabel  = cms.untracked.int32(runOnTTbarMC),
                                     genWeightLabel    = cms.InputTag("genWeight:genWeight"),
                                     genLabel          = cms.InputTag("prunedGenParticles"),
                                     genJetLabel       = cms.InputTag("slimmedGenJets"),
                                     genttbarCatLabel  = cms.InputTag("GenTtbarCategories:genTtbarId"),
                                     genttbarConeCatLabel  = cms.InputTag("catGenTops"),
                                     muonLabel         = cms.InputTag("catMuons"),
                                     electronLabel     = cms.InputTag("catElectrons"),
                                     jetLabel          = cms.InputTag("catJets"),
                                     metLabel         = cms.InputTag("catMETs"),
                                     #metLabel          = cms.InputTag("catMETsNoHF"),
                                     pvLabel           = cms.InputTag("catVertex:nGoodPV"),
                                     puWeightLabel     = cms.InputTag("pileupWeight"),
                                     triggerBits       = cms.InputTag("TriggerResults::HLT"), # Not working yet
                                     triggerObjects    = cms.InputTag("catTrigger"), # Not working yet
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Tree_ttbbLepJets.root')
                                   )


#process.Tracer = cms.Service("Tracer")
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
# process.p = cms.Path(process.pileupWeight*
#                      process.ttbarSingleLepton)
process.p = cms.Path(process.ttbarSingleLepton)
