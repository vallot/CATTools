import FWCore.ParameterSet.Config as cms
process = cms.Process("TopTree")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:catTuple.root'
        '/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v7-6-4_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/160325_073338/0000/catTuple_1.root'
        #'file:/afs/cern.ch/user/b/brochero/CATTuples_August/v7-3-6/cat74/src/CATTools/CatProducer/prod/catTuple-PUPPI-v7-3-6.root' # -- MC PUPPI (v7-3-6)
        #'file:/cms/scratch/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-2_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150805_203807/0000/catTuple_245.root' # -- MC
        #'root://cms-xrdr.sdfarm.kr///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-4_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150810_215031/0000/catTuple_276.root' # -- MC
        #'file:/cms/home/brochero/CATTuples_August/v7-3-4/cat74/src/CATTools/CatAnalyzer/prod/SingleMu-PromptReco_catTuple_44.root'   # -- DATA
        #'root://cms-xrdr.sdfarm.kr///xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-6_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150820_215807/0000/catTuple_108.root'    # -- XROOT test
#        'file:/cms/home/brochero/CATTuples_August/Central-v7-3-6/cat74/src/CATTools/CatAnalyzer/prod/catTuple_108.root' # -- MC
        #'file:/cms/home/brochero/CATTuples_July/cat74/src/CATTools/CatAnalyzer/catTuple_83.root' # -- Data
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_13.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_14.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_16.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_17.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_18.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_19.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_21.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_23.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_24.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_25.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_26.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_27.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_31.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_35.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_36.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_37.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_38.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_41.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_43.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_45.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_6.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_78.root',
       #'/store/user/tjkim/pf/TT_TuneCUETP8M1_13TeV-powheg-pythia8/PAT2CAT/160116_164927/0000/catTuple_9.root',
       #'file:catTuple.root'
    )
)

isData = False

# json file (Only Data)
if isData:
  import FWCore.PythonUtilities.LumiList as LumiList
  process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

process.TopTree = cms.EDAnalyzer('TopAnalyzer',
                                           genTopLabel   = cms.InputTag("catGenTops"),
                                           genLabel      = cms.InputTag("prunedGenParticles"),
                                           muonLabel     = cms.InputTag("catMuons"),
                                           electronLabel = cms.InputTag("catElectrons"),
                                           jetLabel      = cms.InputTag("catJets"),
                                           metLabel      = cms.InputTag("catMETs"),
                                           pvLabel       = cms.InputTag("catVertex:nGoodPV"),
                                           puWeight      = cms.InputTag("pileupWeight"),
                                           genWeightLabel = cms.InputTag("genWeight"),
                                           )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('hep.root')
                                   )


#process.Tracer = cms.Service("Tracer")
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.TopTree)
