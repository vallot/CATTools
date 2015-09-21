#------------------------------------------------------------------
#------------------------------------------------------------------
# Data or MC Sample
runOnMC      = False
# runOnTTbarMC == 0, No ttbar
# runOnTTbarMC == 1, ttbar Signal
# runOnTTbarMC == 2, ttbar Background
runOnTTbarMC = 0
#------------------------------------------------------------------
#------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
process = cms.Process("ttbarSingleLepton")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('ttbarljets')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/user/b/brochero/CATTuples_August/v7-3-6/cat74/src/CATTools/CatProducer/prod/catTuple-PUPPI-v7-3-6.root' # -- MC PUPPI (v7-3-6)
        #'file:/cms/scratch/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150805_203807/0000/catTuple_245.root' # -- MC
        #'root://cms-xrdr.sdfarm.kr///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-4_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150810_215031/0000/catTuple_276.root' # -- MC
        #'file:/cms/home/brochero/CATTuples_August/v7-3-4/cat74/src/CATTools/CatAnalyzer/prod/SingleMu-PromptReco_catTuple_44.root'   # -- DATA
        #'root://cms-xrdr.sdfarm.kr///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150820_215807/0000/catTuple_108.root'    # -- XROOT test
        'file:/cms/home/brochero/CATTuples_August/Central-v7-3-6/cat74/src/CATTools/CatAnalyzer/prod/catTuple_108.root' # -- MC
        #'file:/cms/home/brochero/CATTuples_July/cat74/src/CATTools/CatAnalyzer/catTuple_83.root' # -- Data
    )
)

# json file (Only Data)
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt').getVLuminosityBlockRange()


process.ttbarSingleLepton = cms.EDAnalyzer('TtbarSingleLeptonAnalyzer',
                                           sampleLabel       = cms.untracked.bool(runOnMC),
                                           TTbarSampleLabel  = cms.untracked.int32(runOnTTbarMC),
                                           genLabel      = cms.InputTag("prunedGenParticles"),
                                           muonLabel     = cms.InputTag("catMuons"),
                                           electronLabel = cms.InputTag("catElectrons"),
                                           jetLabel      = cms.InputTag("catJets"),
                                           #metLabel      = cms.InputTag("catMETs"),
                                           metLabel      = cms.InputTag("catMETsNoHF"),
                                           pvLabel       = cms.InputTag("catVertex:nGoodPV"),
                                           puWeight      = cms.InputTag("pileupWeight"),
                                           trigLabel     = cms.InputTag("catTrigger"), # Not working yet
                                           )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('vallot.root')
                                   )


#process.Tracer = cms.Service("Tracer")
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.ttbarSingleLepton)
