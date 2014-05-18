import FWCore.ParameterSet.Config as cms

process = cms.Process("CAT")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Global geometry
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#Data measurements from Summer11
#process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
#process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")

#process.GlobalTag.globaltag = cms.string('START53_V20::All')
# geometry needed for clustering and calo shapes variables
# process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")
# 3 folllowing config files included in RecoEcal.EgammaClusterProducers.geometryForClustering_cff
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

# ES cluster for pi0 discrimination variables
#process.load("RecoEcal.EgammaClusterProducers.preshowerClusterShape_cfi")

# pi0 discrimination variables
#process.load("RecoEcal.EgammaClusterProducers.piZeroDiscriminators_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:PAT.root'))

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.analysis = cms.EDAnalyzer("ObjectProducer",
        myConfig = cms.PSet(
                # Verbosite
                #               0 = muet
                #               1 = Number of evt every 100 evts
                #               2 = Give the functions executed & nof objects build per event
                #               3 = Liste of high level objects (jetss, muons, ...)
                #               4 = List of all  objects 
                #               5 = Debug
                verbosity = cms.untracked.int32(0),

                # used in the electron to see if the magneticfield is taken from DCS or from IDEALMAGFIELDRECORD
                isData = cms.untracked.bool(False),

                # name of output root file
                RootFileName = cms.untracked.string("CAT.root"),

                # What is written to rootuple               
                doHLT = cms.untracked.bool(True),
                doPDFInfo = cms.untracked.bool(True),
                signalGenerator = cms.untracked.string('PYTHIA'),
#               signalGenerator = cms.untracked.string('ALPGEN'),
#               signalGenerator = cms.untracked.string('MADGRAPH'),

                doElectronMC = cms.untracked.bool(True),
                doMuonMC = cms.untracked.bool(True),
                doPhotonMC = cms.untracked.bool(True),
                doJetMC = cms.untracked.bool(True),
                doMETMC = cms.untracked.bool(True),
                doUnstablePartsMC = cms.untracked.bool(True),
                doPrimaryVertex = cms.untracked.bool(True),
                runGeneralTracks = cms.untracked.bool(True),#true only if generalTracks are stored.
                doCaloJet = cms.untracked.bool(False),
                doGenJet = cms.untracked.bool(True),
                doCaloJetId = cms.untracked.bool(False),
                doPFJet = cms.untracked.bool(True),
                doJPTJet = cms.untracked.bool(False),
                doJPTJetId = cms.untracked.bool(False),
                doMuon = cms.untracked.bool(True),
                doElectron = cms.untracked.bool(True),
                doPhoton = cms.untracked.bool(True),
                runSuperCluster = cms.untracked.bool(True),#true only if SuperCluster are stored
                doCaloMET = cms.untracked.bool(False),
                doPFMET = cms.untracked.bool(True),
                doTCMET = cms.untracked.bool(False),
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
                vcaloJetProducer = cms.untracked.vstring("selectedPatJetsAK5Calo"),
                vgenJetProducer = cms.untracked.vstring("ak5GenJetsNoNu"),
                vpfJetProducer = cms.untracked.vstring("selectedPatJetsPF2PAT"),
                vJPTJetProducer = cms.untracked.vstring(""),
                vmuonProducer = cms.untracked.vstring("selectedPatMuonsPF2PAT"),
                velectronProducer = cms.untracked.vstring("selectedPatElectronsPF2PAT"),
                vphotonProducer = cms.untracked.vstring("selectedPatPhotons"),
                CalometProducer = cms.InputTag("patMETs"),
                vpfmetProducer = cms.untracked.vstring("patType1CorrectedPFMetPF2PAT"),
                TCmetProducer = cms.InputTag("patMETsTC"),
                genEventProducer = cms.InputTag("genEvt"),
                generalTrackLabel = cms.InputTag("generalTracks")
        )
)
process.p = cms.Path(process.analysis)

temp = process.dumpPython()
outputfile = file("expanded.py",'w')
outputfile.write(temp)
outputfile.close()
