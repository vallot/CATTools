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
# TTbarCatMC   ==> 0->All ttbar, 1->ttbb, 2->ttcc, 3->ttLF, 4->ttV/H, 5->STFCNC, 6->TTFCNC
options.register('TTbarCatMC', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "TTbarCatMC: 0  default All ttbar events")
options.parseArguments()

print "User JSON file: " + str(options.UserJSON)
print "runOnTTbarMC: "   + str(options.runOnTTbarMC)
print "TTbarCatMC: "     + str(options.TTbarCatMC)
#print "PU Map: "         + str(options.PUMap)
#------------------------------------------------------------------
#------------------------------------------------------------------

process = cms.Process("fcncLepJets")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('fcncLepJets')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit = cms.untracked.int32(-1)
#)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        'root://cluster142.knu.ac.kr:1094///store/group/CAT/V10_3/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/V10_3_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/200725_030313/0000/catTuple_1.root'
        #'root://cluster142.knu.ac.kr:1094//store/group/CAT/V10_3/ST_FCNC-TH_Tleptonic_HTobb_eta_hut-MadGraph5-pythia8/V10_3_RunIIAutumn18MiniAOD-tauDecays_102X_upgrade2018_realistic_v15-v1/201209_124929/0000/catTuple_1.root'
        #'root://cluster142.knu.ac.kr:1094//store/group/CAT/V10_3/TT_FCNC-TtoHJ_aTleptonic_HTobb_eta_hut_TuneCP5-MadGraph5-pythia8/V10_3_RunIIAutumn18MiniAOD-tauDecays_102X_upgrade2018_realistic_v15-v1/191230_005512/0000/catTuple_1.root'
        #'root://cluster142.knu.ac.kr:1094///store/group/CAT/V10_3/SingleMuon/V10_3_Run2018B-17Sep2018-v1/191228_195656/0000/catTuple_2.root'
        )
)

# PUReWeight
process.load("CATTools.CatProducer.pileupWeight_cff")
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.pileupWeight.weightingMethod = "RedoWeight"
#process.pileupWeight.pileupMC = pileupWeightMap[options.PUMap]
process.pileupWeight.pileupMC = pileupWeightMap["2018_25ns_MC"]
process.pileupWeight.pileupRD = pileupWeightMap["Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON"]
process.pileupWeight.pileupUp = pileupWeightMap["Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_Up"]
process.pileupWeight.pileupDn = pileupWeightMap["Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_Dn"]

# json file (Only Data)
if options.UserJSON:
    # /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco
    print "Running data.... Including JSON File."
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '../../CatProducer/data/LumiMask/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt').getVLuminosityBlockRange()

#Load MET filters for MC
process.load("CATTools.CatAnalyzer.filters_cff")
#process.filterRECOMC.doFilter = True #double check
# Lepton Scale Factors
from CATTools.CatAnalyzer.leptonSF_cff import *
# GEN Weights
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
process.flatGenWeights.keepFirstOnly = False

process.fcncLepJets = cms.EDAnalyzer('fcncLepJetsAnalyzer',
                                     TTbarSampleLabel  = cms.untracked.int32(options.runOnTTbarMC),
                                     TTbarCatLabel     = cms.untracked.int32(options.TTbarCatMC),
                                     trigMuFilters     = cms.InputTag("filterTrigMU"),
                                     trigElFilters     = cms.InputTag("filterTrigEL"),
                                     trigElHTFilters   = cms.InputTag("filterTrigELHT"),
                                     recoFiltersMC     = cms.InputTag("filterRECOMC"),
                                     recoFilters       = cms.InputTag("filterRECO"),
                                     # JES file
                                     JESUncFile        = cms.untracked.string("RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt"),
                                     # Input Tags
                                     genWeightLabel    = cms.InputTag("flatGenWeights"),
                                     allweightLabel    = cms.InputTag("genWeight"),
                                     genLabel          = cms.InputTag("prunedGenParticles"),
                                     genJetLabel       = cms.InputTag("slimmedGenJets"),
                                     genHiggsCatLabel  = cms.InputTag("GenTtbarCategories:genTtbarId"),
                                     genttbarCatLabel  = cms.InputTag("catGenTops"),
                                     muonLabel         = cms.InputTag("catMuons"),
                                     muonIdSF          = muonSFTightId102X,
                                     muonIsoSF         = muonSFTightIso102X,
                                     muonTrgSF         = trigSF_IsoMu24,
                                     electronLabel     = cms.InputTag("catElectrons"),
                                     elecIdSF          = electronSFCutBasedTight102X,
                                     elecRecoSF        = electronSFReco102X,
                                     elecZvtxSF        = electronSFHLTZvtx102X,
                                     elecTrgSF         = trigSF_El32_El28HT150_ttH_legacy18_v1,
                                     jetLabel          = cms.InputTag("catJets"),
                                     metLabel          = cms.InputTag("catMETs"),
                                     pvLabel           = cms.InputTag("catVertex:nGoodPV"),
                                     puWeightLabel     = cms.InputTag("pileupWeight"),
                                     triggerBits       = cms.InputTag("TriggerResults"), 
                                     triggerObjects    = cms.InputTag("catTrigger"), 
                                     JetMother         = cms.InputTag("genJetHadronFlavour:ancestors"),
                                     nTrueVertLabel    = cms.InputTag("pileupWeight:nTrueInteraction")
                                     )

process.fcncLepJetsQCD = process.fcncLepJets.clone(
    doLooseLepton = cms.untracked.bool(True),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Tree_fcncLepJets.root')
                                   )

process.p = cms.Path(process.filterRECOMC + process.filterRECO +
                     process.filterTrigMU + process.filterTrigEL + process.filterTrigELJET + process.filterTrigELHT +
                     process.flatGenWeights +
                     process.pileupWeight +
                     process.fcncLepJets)#+ process.fcncLepJetsQCD)
