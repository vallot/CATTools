from PhysicsTools.PatAlgos.patTemplate_cfg import *
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

#process.source.fileNames = filesRelValProdTTbarAODSIM
# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_0_0/RelValTTbar_13/GEN-SIM-RECO/PU50ns_POSTLS170_V4-v2/00000/265B9219-FF98-E311-BF4A-02163E00EA95.root')
)

process.maxEvents.input = 100


process.analysis = cms.EDAnalyzer("ObjectProducer",
        myConfig = cms.PSet(
                # Verbosite
                #               0 = muet
                #               1 = Number of evt every 100 evts
                #               2 = Give the functions executed & nof objects build per event
                #               3 = Liste of high level objects (jetss, muons, ...)
                #               4 = List of all  objects 
                #               5 = Debug
                verbosity = cms.untracked.int32(1),

                # used in the electron to see if the magneticfield is taken from DCS or from IDEALMAGFIELDRECORD
                isData = cms.untracked.bool(False),

                # name of output root file
                RootFileName = cms.untracked.string("CAT.root"),

                # What is written to rootuple               
                doHLT = cms.untracked.bool(False),
                doPDFInfo = cms.untracked.bool(False),
                signalGenerator = cms.untracked.string('PYTHIA'),
#               signalGenerator = cms.untracked.string('ALPGEN'),
#               signalGenerator = cms.untracked.string('MADGRAPH'),

                doElectronMC = cms.untracked.bool(False),
                doMuonMC = cms.untracked.bool(False),
                doPhotonMC = cms.untracked.bool(False),
                doJetMC = cms.untracked.bool(False),
                doMETMC = cms.untracked.bool(False),
                doUnstablePartsMC = cms.untracked.bool(False),
                doPrimaryVertex = cms.untracked.bool(True),
                runGeneralTracks = cms.untracked.bool(False),#true only if generalTracks are stored.
                doGenJet = cms.untracked.bool(False),
                doPFJet = cms.untracked.bool(True),
                doMuon = cms.untracked.bool(True),
                doElectron = cms.untracked.bool(True),
                doPhoton = cms.untracked.bool(True),
                runSuperCluster = cms.untracked.bool(True),#true only if SuperCluster are stored
                doPFMET = cms.untracked.bool(True),
                doGenEvent = cms.untracked.bool(False),#put on False when running non-ttbar or when running toptree from reco
                doNPGenEvent = cms.untracked.bool(False),#put on True when running New Physics sample
                doSpinCorrGen = cms.untracked.bool(False),#put on True only if you need SpinCorrelation Variables
                doSemiLepEvent = cms.untracked.bool(False),#put on True only if you need TtSemiLeptonicEvent Collection exist in PAT-uples (L2)

                conversionLikelihoodWeightsFile = cms.untracked.string('RecoEgamma/EgammaTools/data/TMVAnalysis_Likelihood.weights.txt'),

                # Draw MC particle tree
                drawMCTree = cms.untracked.bool(False),
                mcTreePrintP4 = cms.untracked.bool(False),
                mcTreePrintPtEtaPhi = cms.untracked.bool(False),
                mcTreePrintVertex = cms.untracked.bool(False),
                mcTreePrintStatus = cms.untracked.bool(False),
                mcTreePrintIndex = cms.untracked.bool(False),
                mcTreeStatus = cms.untracked.vint32( 3 ),       # accepted status codes


                # MC particles acceptance cuts
                electronMC_etaMax = cms.double(3.0),
                electronMC_ptMin = cms.double(2.0),
                muonMC_etaMax = cms.double(3.0),
                muonMC_ptMin = cms.double(0.0),
                jetMC_etaMax = cms.double(6.0),
                jetMC_ptMin = cms.double(5.0),
                useEventCounter = cms.untracked.bool( True ),
                filters = cms.untracked.vstring(
                                             'prePathCounter',
                                             'postPathCounter'
                                               )
        ),

        producersNames = cms.PSet(
                hltProducer1st = cms.InputTag("TriggerResults","","HLT"),
                hltProducer2nd = cms.InputTag("TriggerResults","","RECO"),
                hltProducer3rd = cms.InputTag("TriggerResults","","REDIGI"),
                hltProducer4th = cms.InputTag("TriggerResults","","REDIGI311X"),
                pileUpProducer = cms.InputTag("addPileupInfo"),
                genParticlesProducer = cms.InputTag("genParticles"),
                primaryVertexProducer = cms.InputTag("goodOfflinePrimaryVertices"),
                vgenJetProducer = cms.untracked.vstring("ak5GenJetsNoNu"),
                vpfJetProducer = cms.untracked.vstring("selectedPatJets"),
                vmuonProducer = cms.untracked.vstring("selectedPatMuons"),
                velectronProducer = cms.untracked.vstring("selectedPatElectrons"),
                vphotonProducer = cms.untracked.vstring("selectedPatPhotons"),
                vpfmetProducer = cms.untracked.vstring("patMETs"),
                genEventProducer = cms.InputTag("genEvt"),
                generalTrackLabel = cms.InputTag("generalTracks")
        )
)

process.p = cms.Path(process.analysis)
