import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
#'/store/user/youn/cat710_phy14_ttbar_2025_aod/catTuple_972.root'
#'file:catTuple_1.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_1.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_2.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_3.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_4.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_5.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_6.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_7.root',
#'/store/group/CAT/MuonEG/v7-3-4_Run2015B-PromptReco-v1/150810_214338/0000/catTuple_8.root',

#'/store/user/youn/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20150722_231402/150722_141433/0000/catTuple_1.root'

#'/store/group/CAT/MuonEG/desySync_Run2015B-PromptReco-v1/150905_193922/0000/catTuple_12.root',
#'/store/group/CAT/MuonEG/desySync_Run2015B-17Jul2015-v1/150905_193954/0000/catTuple_1.root',
#kisti:
#'file:/cms/scratch/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_1.root',
#'file:/cms/scratch/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_2.root',
#'file:/cms/scratch/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_3.root',
#'file:/cms/scratch/geonmo/cattuple/catTuple_1.root',
#'file:/cms/scratch/geonmo/cattuple/catTuple_2.root',

#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_99.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_28.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_81.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_44.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_87.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_59.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_72.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_86.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_90.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_26.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_47.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_55.root',
#'/store/group/CAT/DoubleEG/v7-3-6_Run2015B-PromptReco-v1/150922_133632/0000/catTuple_45.root',

#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_73.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_34.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_13.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_67.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_54.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_33.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_22.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_47.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_72.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_56.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_39.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_11.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_69.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_15.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_5.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_83.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_77.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_37.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_65.root',
#'/store/group/CAT/MuonEG/v7-4-1_Run2015D-PromptReco-v3/151002_181515/0000/catTuple_41.root',

'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_52.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_71.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_16.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_17.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_30.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_31.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_27.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_39.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_47.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_55.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_59.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_15.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_22.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_10.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_20.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_28.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_26.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_42.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_44.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_25.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_49.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_57.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_48.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_77.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_45.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_64.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_72.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_76.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_58.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_68.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_14.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_19.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_12.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_13.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_11.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_34.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_38.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_23.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_46.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_53.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_54.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_36.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_63.root',
'/store/group/CAT/MuonEG/v7-4-3_Run2015C-PromptReco-v1/151009_212654/0000/catTuple_60.root',

#uos : gate2
#'/cms/store/group/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_1.root',
#'/cms/store/group/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_2.root',
#'/cms/store/group/CAT/MuonEG/v7-3-0_Run2015B-PromptReco-v1/150720_060935/0000/catTuple_3.root',

#'file:/cms/home/youn/work/cattool/tag711/cat/src/CATTools/CatProducer/prod/catTuple.root'
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_1.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_2.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_3.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_4.root',
#          '/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/141219_091640/0000/catTuple_5.root',
      )
)

from FWCore.PythonUtilities import LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'desy.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2_1.txt').getVLuminosityBlockRange()#jul
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2_2.txt').getVLuminosityBlockRange()#prompt
process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-257599_13TeV_PromptReco_Collisions15_25ns_JSON_C.txt').getVLuminosityBlockRange()

#process.source.lumisToProcess = LumiList.LumiList(filename = 'json_DCSONLY_Run2015B.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON_MuonPhys.txt').getVLuminosityBlockRange()

process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring("nEventsTotal"), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
                                
#    bool = cms.PSet(
#        CSCTightHaloFilter = cms.PSet(src = cms.InputTag("catTrigger","CSCTightHaloFilter")),
#        EcalDeadCellTriggerPrimitiveFilter =cms.PSet(src = cms.InputTag("catTrigger","EcalDeadCellTriggerPrimitiveFilter")),
#        HBHENoiseFilter = cms.PSet(src = cms.InputTag("catTrigger","HBHENoiseFilter")),
#        eeBadScFilter = cms.PSet(src = cms.InputTag("catTrigger","eeBadScFilter")),
#        goodVertices = cms.PSet(src = cms.InputTag("catTrigger","goodVertices")),
#    ),
                                
    int = cms.PSet(
        modes = cms.PSet(src = cms.InputTag("partonTop","modes")),
        channel = cms.PSet(src = cms.InputTag("partonTop","channel")),
        
        #HLTDoubleEle33CaloIdLGsfTrkIdVL = cms.PSet(src = cms.InputTag("catTrigger","HLTDoubleEle33CaloIdLGsfTrkIdVL")),
        #HLTEle12CaloIdLTrackIdLIsoVL = cms.PSet(src = cms.InputTag("catTrigger","HLTEle12CaloIdLTrackIdLIsoVL")),
        #HLTEle16Ele12Ele8CaloIdLTrackIdL = cms.PSet(src = cms.InputTag("catTrigger","HLTEle16Ele12Ele8CaloIdLTrackIdL")),
        #HLTEle17CaloIdLTrackIdLIsoVL = cms.PSet(src = cms.InputTag("catTrigger","HLTEle17CaloIdLTrackIdLIsoVL")),
        #HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ = cms.PSet(src = cms.InputTag("catTrigger","HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ")),
        #HLTEle23Ele12CaloIdLTrackIdLIsoVL = cms.PSet(src = cms.InputTag("catTrigger","HLTEle23Ele12CaloIdLTrackIdLIsoVL")),
        #HLTEle23Ele12CaloIdLTrackIdLIsoVLDZ = cms.PSet(src = cms.InputTag("catTrigger","HLTEle23Ele12CaloIdLTrackIdLIsoVLDZ")),
        #HLTEle27eta2p1WPLooseGsfTriCentralPFJet30 = cms.PSet(src = cms.InputTag("catTrigger","HLTEle27eta2p1WPLooseGsfTriCentralPFJet30")),
        #HLTMu17Mu8DZ = cms.PSet(src = cms.InputTag("catTrigger","HLTMu17Mu8DZ")),
        #HLTMu17TkMu8DZ = cms.PSet(src = cms.InputTag("catTrigger","HLTMu17TkMu8DZ")),
        #HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL = cms.PSet(src = cms.InputTag("catTrigger","HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL")),
        #HLTMu17TrkIsoVVLMu8TrkIsoVVL = cms.PSet(src = cms.InputTag("catTrigger","HLTMu17TrkIsoVVLMu8TrkIsoVVL")),
        #HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ = cms.PSet(src = cms.InputTag("catTrigger","HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ")),
        #HLTMu17TrkIsoVVLTkMu8TrkIsoVVL = cms.PSet(src = cms.InputTag("catTrigger","HLTMu17TrkIsoVVLTkMu8TrkIsoVVL")),
        #HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL = cms.PSet(src = cms.InputTag("catTrigger","HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL")),
                   
        nGoodPV = cms.PSet(src = cms.InputTag("catVertex","nGoodPV")),
        nPV = cms.PSet(src = cms.InputTag("catVertex","nPV")),
           
                   
                   
        #nVertex   = cms.PSet(src = cms.InputTag("recoEventInfo","pvN")),
        #HLTDoubleMu = cms.PSet(src = cms.InputTag("recoEventInfo","HLTDoubleMu")),
        #HLTDoubleEl = cms.PSet(src = cms.InputTag("recoEventInfo","HLTDoubleEl")),
        #HLTMuEl = cms.PSet(src = cms.InputTag("recoEventInfo","HLTMuEl")),
        #HLTSingleMu = cms.PSet(src = cms.InputTag("recoEventInfo","HLTSingleMu")),
        #HLTSingleEl = cms.PSet(src = cms.InputTag("recoEventInfo","HLTSingleEl")),
    ),
    
    doubles = cms.PSet(
        pdfWeight = cms.PSet(src = cms.InputTag("pdfWeight")),
    ),
    cands = cms.PSet(
        muon = cms.PSet(
            src = cms.InputTag("catMuons"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                e   = cms.string("energy"),
                #relIso = cms.string("relIso"),
                relIso03 = cms.string("relIso(0.3)"),
                relIso04 = cms.string("relIso(0.4)"),
                isTracker = cms.string("isTrackerMuon"),
                isGlobal = cms.string("isGlobalMuon"),
                isLoose = cms.string("isLooseMuon"),
                isTight = cms.string("isTightMuon"),
                dxy = cms.string("dxy"),
                dz = cms.string("dz"),
                q = cms.string("charge"),
                #matched = cms.string("mcMatched"),
            ),
            selections = cms.untracked.PSet(),
        ),
        electrons = cms.PSet(
            src = cms.InputTag("catElectrons"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                e   = cms.string("energy"),
                idMedium = cms.string("electronID('cutBasedElectronID-Spring15-25ns-V1-standalone-medium')"),
                idVeto = cms.string("electronID('cutBasedElectronID-Spring15-25ns-V1-standalone-veto')"),
                #mva = cms.string("electronID('mvaTrigV0')"),
                relIso03 = cms.string("relIso(0.3)"),
                relIso04 = cms.string("relIso(0.4)"),
                scEta = cms.string("scEta"),
                q = cms.string("charge"),
                passConversionVeto = cms.string("passConversionVeto"),
                isPF = cms.string("isPF"),
                               
            ),
            selections = cms.untracked.PSet(
                                            #isPassBaseId = cms.string("passConversionVeto && isPF && gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS') <= 0"),
                #isPassBaseId = cms.string("passConversionVeto && isPF"),
            ),
        ),
        jets = cms.PSet(
            src = cms.InputTag("catJets"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                e   = cms.string("energy"),
                vtxMass = cms.string("vtxMass"),
                CSVInclV2 = cms.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),
                #CSVInclV2 = cms.string("bDiscriminator('combinedInclusiveSecondaryVertexV2BJetTags')"),
                partonFlavour = cms.string("partonFlavour"),
                hadronFlavour = cms.string("hadronFlavour"),
                isLoose = cms.string("LooseId"),
                isPFId = cms.string("pileupJetId"),
            ),
            selections = cms.untracked.PSet(
                #isLoose = cms.string("LooseId"),
                #isPFId = cms.string("pileupJetId"),
            ),
        ),
                     
        jetsPuppi = cms.PSet(
            src = cms.InputTag("catJetsPuppi"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                e   = cms.string("energy"),
                vtxMass = cms.string("vtxMass"),
                CSVInclV2 = cms.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),
                #partonFlavour = cms.string("partonFlavour"),
                #hadronFlavour = cms.string("hadronFlavour"),
                isLoose = cms.string("LooseId"),
                isPFId = cms.string("pileupJetId"),
            ),
            selections = cms.untracked.PSet(
#isLoose = cms.string("LooseId"),
                #isPFId = cms.string("pileupJetId"),
            ),
        ),
             
        met = cms.PSet(
            src = cms.InputTag("catMETs"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
        noHFmet = cms.PSet(
            src = cms.InputTag("catMETsNoHF"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),

        metPfMva = cms.PSet(
            src = cms.InputTag("catMETsPfMva"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
                     
        metPuppi = cms.PSet(
            src = cms.InputTag("catMETsPuppi"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                phi = cms.string("phi"),
            ),
            selections = cms.untracked.PSet(),
        ),
             
                     
        slimmedGenJets = cms.PSet(
            src = cms.InputTag("slimmedGenJets",""),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                #pdgId = cms.string("pdgId"),
                #q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        partonTop = cms.PSet(
            src = cms.InputTag("partonTop"),
            #modes = cms.string("modes"),
            #channel = cms.string("channel"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTopJet = cms.PSet(
            src = cms.InputTag("pseudoTop","jets"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTopLepton = cms.PSet(
            src = cms.InputTag("pseudoTop","leptons"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTopNu = cms.PSet(
            src = cms.InputTag("pseudoTop","neutrinos"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                #m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                #q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        pseudoTop = cms.PSet(
            src = cms.InputTag("pseudoTop"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

#process.load("CATTools.CatProducer.pseudoTop_cff")
#process.p = cms.Path(
#    process.nEventsTotal*
#    process.partonTop*
#    process.ntuple
#)


process.p = cms.Path(
    process.nEventsTotal*
    process.ntuple
)

