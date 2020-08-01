import FWCore.ParameterSet.Config as cms
#------------------------------------------------------------------
#------------------------------------------------------------------
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# JSON
options.register('UserJSON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "UserJSON: Fault  default")
# runOnTTbarMC ==> 0->No ttbar, 1->ttbar Signal
options.register('runOnTTbarMC', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "runOnTTbarMC: 0  default No ttbar sample")
options.parseArguments()

print "User JSON file: " + str(options.UserJSON)
print "runOnTTbarMC: "   + str(options.runOnTTbarMC)
#print "PU Map: "         + str(options.PUMap)
#------------------------------------------------------------------
#------------------------------------------------------------------

process = cms.Process("test")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        'file:/pnfs/knu.ac.kr/data/cms/store/group/CAT/V8_1/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/V8_1_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/200725_154813/0000/catTuple_1.root'
        #'file:/pnfs/knu.ac.kr/data/cms/store/group/CAT/V8_1/SingleMuon/V8_1_Run2016B-17Jul2018_ver2-v1/200414_062856/0000/catTuple_1.root'
        )
)

# json file (Only Data)
if options.UserJSON:
    # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco
    print "Running data.... Including JSON File."
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '../../CatProducer/data/LumiMask/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt').getVLuminosityBlockRange()

#Load MET filters for MC
process.load("CATTools.CatAnalyzer.filters_cff")
#process.filterRECOMC.doFilter = True #double check
# Lepton Scale Factors
from CATTools.CatAnalyzer.leptonSF_cff import *
# GEN Weights
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
process.flatGenWeights.keepFirstOnly = False

process.test = cms.EDAnalyzer('testAnalyzer',
                                     TTbarSampleLabel  = cms.untracked.int32(options.runOnTTbarMC),
                                     trigMuFilters     = cms.InputTag("filterTrigMU"),
                                     trigElFilters     = cms.InputTag("filterTrigEL"),
                                     trigElHTFilters   = cms.InputTag("filterTrigELHT"),
                                     recoFiltersMC     = cms.InputTag("filterRECOMC"),
                                     recoFilters       = cms.InputTag("filterRECO"),
                                     # Input Tags
                                     genWeightLabel    = cms.InputTag("flatGenWeights"),
                                     genLabel          = cms.InputTag("prunedGenParticles"),
                                     genJetLabel       = cms.InputTag("slimmedGenJets"),
                                     genHiggsCatLabel  = cms.InputTag("GenTtbarCategories:genTtbarId"),
                                     genttbarCatLabel  = cms.InputTag("catGenTops"),
                                     muonLabel         = cms.InputTag("catMuons"),
                                     electronLabel     = cms.InputTag("catElectrons"),
                                     jetLabel          = cms.InputTag("catJets"),
                                     metLabel          = cms.InputTag("catMETs"),
                                     pvLabel           = cms.InputTag("catVertex:nGoodPV"),
                                     puWeightLabel     = cms.InputTag("pileupWeight"),
                                     triggerBits       = cms.InputTag("TriggerResults"), 
                                     triggerObjects    = cms.InputTag("catTrigger"), 
                                     JetMother         = cms.InputTag("genJetHadronFlavour:ancestors"),
                                     nTrueVertLabel    = cms.InputTag("pileupWeight:nTrueInteraction")
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Tree_test.root')
                                   )

process.p = cms.Path(process.filterRECOMC + process.filterRECO +
                     process.filterTrigMU + process.filterTrigEL +
                     process.flatGenWeights +
                     process.test)
