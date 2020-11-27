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
        'root://cluster142.knu.ac.kr:1094//store/group/CAT/V9_7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/V9_7_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/200723_184504/0000/catTuple_1.root'
        #'root://cluster142.knu.ac.kr:1094//store/group/CAT/V9_7/SingleMuon/V9_7_Run2017B-31Mar2018-v1/200723_112209/0000/catTuple_1.root'
        )
)

# json file (Only Data)
if options.UserJSON:
    # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco
    print "Running data.... Including JSON File."
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '../../CatProducer/data/LumiMask/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt').getVLuminosityBlockRange()

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
                     process.filterTrigMU + process.filterTrigEL + process.filterTrigELJET + process.filterTrigELHT +
                     process.flatGenWeights +
                     process.test)
