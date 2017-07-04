import FWCore.ParameterSet.Config as cms
#------------------------------------------------------------------
#------------------------------------------------------------------
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# JSON
options.register('UserJSON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "UserJSON: Fault  default")
# runOnTTbarMC ==> 0->No ttbar, 1->ttbar Signal, 2->ttbar Background
options.register('runOnTTbarMC', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "runOnTTbarMC: 0  default No ttbar sample")
# TTbarCatMC   ==> 0->All ttbar, 1->ttbb, 2->ttbj, 3->ttcc, 4->ttLF, 5->tt, 6->ttjj
options.register('TTbarCatMC', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "TTbarCatMC: 0  default All ttbar events")
options.parseArguments()

print "User JSON file: " + str(options.UserJSON)
print "runOnTTbarMC: "   + str(options.runOnTTbarMC)
print "TTbarCatMC: "     + str(options.TTbarCatMC)
#------------------------------------------------------------------
#------------------------------------------------------------------

process = cms.Process("ttbbLepJets")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
# process.MessageLogger.cerr.threshold = 'INFO'
# process.MessageLogger.categories.append('ttbbLepJets')
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)
# )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",

     # fileNames = cms.untracked.vstring()
     fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/user/b/brochero/brochero_WorkArea/CATTuples_2017/CMSSW_8_0_20/src/CATTools/CatAnalyzer/prod/catTuple_ttbarLepJetsPowhegv803.root',
        #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-4_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170113_045423/0000/catTuple_1.root',
        'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_104856/0000/catTuple_106.root',
        #'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/SingleElectron/v8-0-3_Run2016E-23Sep2016-v1/161207_134721/0000/catTuple_1.root'
        )
)
# from CATTools.Validation.commonTestInput_cff import commonTestCATTuples
# process.source.fileNames = commonTestCATTuples["sig"]

# PUReWeight
process.load("CATTools.CatProducer.pileupWeight_cff")
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.pileupWeight.weightingMethod = "RedoWeight"
process.pileupWeight.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.pileupWeight.pileupRD = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON"]
process.pileupWeight.pileupUp = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Up"]
process.pileupWeight.pileupDn = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Dn"]

# json file (Only Data)
if options.UserJSON:
    # ReReco JSON file taken from: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt
    print "Running data.... Including JSON File."
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

# Lepton Scale Factors
from CATTools.CatAnalyzer.leptonSF_cff import *
# GEN Weights
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
# CSV Scale Factors
process.load("CATTools.CatAnalyzer.csvWeights_cfi")
process.csvWeights.minPt = 30   # Same cuts than jet selection
process.csvWeights.maxEta = 2.4

process.ttbbLepJets = cms.EDAnalyzer('ttbbLepJetsAnalyzer',
                                     TTbarSampleLabel  = cms.untracked.int32(options.runOnTTbarMC),
                                     TTbarCatLabel     = cms.untracked.int32(options.TTbarCatMC),
                                     # Skim: Cut in number of (gen) jets 
                                     Skim_N_Jets       = cms.untracked.uint32(0),
                                     # Constrain in Kin. Fitter using CSV position
                                     KFUsebtag         = cms.untracked.bool(True),
                                     CSVPosConKF       = cms.untracked.bool(True),
                                     # TriggerNames
                                     triggerNameDataEl = cms.untracked.vstring("HLT_Ele27_eta2p1_WPTight_Gsf_v","HLT_Ele32_eta2p1_WPTight_Gsf_v"), 
                                     triggerNameDataMu = cms.untracked.vstring("HLT_IsoMu24_v","HLT_IsoTkMu24_v"), 
                                     triggerNameMCEl   = cms.untracked.vstring("HLT_Ele32_eta2p1_WPTight_Gsf_v8"), 
                                     triggerNameMCMu   = cms.untracked.vstring("HLT_IsoMu24_v4","HLT_IsoTkMu24_v4"), 
                                     # Input Tags
                                     genWeightLabel    = cms.InputTag("flatGenWeights"),
                                     genLabel          = cms.InputTag("prunedGenParticles"),
                                     genJetLabel       = cms.InputTag("slimmedGenJets"),
                                     csvWeightLabel    = cms.InputTag("csvWeights"),
                                     genHiggsCatLabel  = cms.InputTag("GenTtbarCategories:genTtbarId"),
                                     genttbarCatLabel  = cms.InputTag("catGenTops"),
                                     muonLabel         = cms.InputTag("catMuons"),
                                     muonSF            = muonSFTight,
                                     electronLabel     = cms.InputTag("catElectrons"),
                                     elecSF            = electronSFCutBasedIDMediumWP,
                                     jetLabel          = cms.InputTag("catJets"),
                                     metLabel          = cms.InputTag("catMETs"),
                                     pvLabel           = cms.InputTag("catVertex:nGoodPV"),
                                     puWeightLabel     = cms.InputTag("pileupWeight"),
                                     triggerBits       = cms.InputTag("TriggerResults"), 
                                     triggerObjects    = cms.InputTag("catTrigger"), 
                                     JetMother         = cms.InputTag("genJetHadronFlavour:ancestors"),
                                     )

process.ttbbLepJetsQCD = process.ttbbLepJets.clone(
    doLooseLepton = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Tree_ttbbLepJets.root')
                                   )


# process.Tracer = cms.Service("Tracer")
# process.dump=cms.EDAnalyzer('EventContentAnalyzer')
# process.p = cms.Path(process.demo*process.dump)
# process.p = cms.Path(process.pileupWeight*
#                      process.ttbarSingleLepton)
process.p = cms.Path(process.flatGenWeights +
                     process.csvWeights +
                     process.pileupWeight +
                     process.ttbbLepJets + process.ttbbLepJetsQCD)
