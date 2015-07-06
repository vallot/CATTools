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
        #'file:/cms/home/brochero/CMSSW_7_3_1/src/CATTools/CommonTools/tag711_phy14_2025_aod/DY_catTuple_078.root',
        #'file:/cms/home/brochero/CMSSW_7_3_1/src/CATTools/CommonTools/tag711_phy14_2025_aod/ttbar_catTuple_535.root',
        #'file:/cms/home/brochero/CMSSW_7_3_1/src/CATTools/CommonTools/tag711_phy14_2025_aod/wj_catTuple_032.root',
        #'file:/cms/home/brochero/CMSSW_7_3_1/src/CATTools/CommonTools/tag711_phy14_2025_aod/tw_catTuple_032.root',
        #'root://cms-xrdr.sdfarm.kr//cms/data/xrd/store/user/youn/nt_tag711_ttbar_2025_aod_v2/ntuple_917.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_000.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_001.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_002.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_003.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_004.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_005.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_006.root',
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_007.root', 
        # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_008.root',
         #'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/youn/tag710_phy14_wj_2025_aod/catTuple_009.root',
        #'root://cms-xrdr.sdfarm.kr///cms/home/brochero/CMSSW_7_3_1/src/CATTools/CommonTools/crab_CAT_QCD_MuEnriched_MiniAOD/results/catTuple_34.root',
       # 'file:/cms/home/brochero/CMSSW_7_3_1/src/CATTools/CommonTools/crab_CAT_QCD_MuEnriched_MiniAOD/results/catTuple_34.root'
       # 'root://cms-xrdr.sdfarm.kr///cms/data/xrd/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/cat74v1_Phys14DR-PU20bx25_PHYS14_25_V1-v1/150616_191303/0000/catTuple_96.root',
        'file:/cms/data/xrd/store/user/jlee/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/cat74v1_Phys14DR-PU20bx25_PHYS14_25_V1-v1/150616_191303/0000/catTuple_96.root',
       #'file:../../CatProducer/prod/catTuple.root' # this is for test using catTuple.root produced locally.
       # '',
    )
)

process.ttbarSingleLepton = cms.EDAnalyzer('TtbarSingleLeptonAnalyzer',
           muonLabel = cms.InputTag("catMuons"),
           electronLabel = cms.InputTag("catElectrons"),
           jetLabel = cms.InputTag("catJets"),
           metLabel = cms.InputTag("catMETs"),
           pvLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
           puWeight = cms.InputTag("pileupWeight")
)

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('vallot.root')
                                   )


#process.Tracer = cms.Service("Tracer")
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.ttbarSingleLepton)
