# for all pflow stuff, we use PAT+PF2PAT now. put the following in your PAT cfg to have the PF2PAT running

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences. It is currently
# not possible to run PF2PAT+PAT and standart PAT at the same time

#from PhysicsTools.PatAlgos.tools.pfTools import *

# An empty postfix means that only PF2PAT is run,  
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
## collections have standard names + postfix (e.g. patElectronPFlow)

#postfix = "PF"
#usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix)
#getattr(process, "patElectrons"+postfix).embedGenMatch = False
#getattr(process, "patMuons"+postfix).embedGenMatch = False

# finally you need to add getattr(process,"patPF2PATSequence"+postfix) to your path.

import FWCore.ParameterSet.Config as cms

process = cms.Process("NewProcess")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Global geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START311_V2::All')
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

process.source = cms.Source("PoolSource",
#        duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
	fileNames = cms.untracked.vstring('file:test311X_PAT_genEvent.root')
	#fileNames = cms.untracked.vstring('/store/data/CRAFT09/Cosmics/RAW-RECO/SuperPointing-CRAFT09_R_V4_CosmicsSeq_v1/0009/763782DB-DCB9-DE11-A238-003048678B30.root')
	#fileNames = cms.untracked.vstring('file:/user/echabert/CMSSW/CMSSW_2_2_3/src/TopQuarkAnalysis/TopEventProducers/test/toto2.root')
#	fileNames = cms.untracked.vstring('dcap:///pnfs/iihe/cms/store/user/blyweert/temp/CRAFT09_RERECO/CRAFT09_SuperPointing_RERECO_CMSSW330_1.root')
)

process.analysis = cms.EDAnalyzer("CatProducer",
	myConfig = cms.PSet(
		# Verbosite
		# 		0 = muet
 		# 		1 = Number of evt every 100 evts
 		# 		2 = Give the functions executed & nof objects build per event
 		# 		3 = Liste of high level objects (jetss, muons, ...)
 		# 		4 = List of all  objects 
		# 		5 = Debug
 		verbosity = cms.untracked.int32(0),

		# used in the electron to see if the magneticfield is taken from DCS or from IDEALMAGFIELDRECORD
		# also used in the JetAnalyzer to see if L2L3Residual correction needs to be stored
		isData = cms.untracked.bool(False),

		# name of output root file
		RootFileName = cms.untracked.string('test311X_TOPTREE.root'),

		# What is written to rootuple		    
		doHLT = cms.untracked.bool(False),
		doMC = cms.untracked.bool(False),
		doPDFInfo = cms.untracked.bool(True),
		signalGenerator = cms.untracked.string('PYTHIA'),
#		signalGenerator = cms.untracked.string('ALPGEN'),
#		signalGenerator = cms.untracked.string('MADGRAPH'),

		doElectronMC = cms.untracked.bool(False),
		doMuonMC = cms.untracked.bool(False),
		doJetMC = cms.untracked.bool(False),
		doMETMC = cms.untracked.bool(False),
		doUnstablePartsMC = cms.untracked.bool(False),
		doPrimaryVertex = cms.untracked.bool(False),
		runGeneralTracks = cms.untracked.bool(False),#true only if generalTracks are stored.
		doCaloJet = cms.untracked.bool(True),
		doGenJet = cms.untracked.bool(True),
		doCaloJetId = cms.untracked.bool(True),
		doPFJet = cms.untracked.bool(True),
    doJPTJet = cms.untracked.bool(False),
		doJPTJetId = cms.untracked.bool(True),
		doMuon = cms.untracked.bool(True),
		doElectron = cms.untracked.bool(True),
		runSuperCluster = cms.untracked.bool(False),#true only if SuperCluster are stored
		doCaloMET = cms.untracked.bool(True),
		doPFMET = cms.untracked.bool(False),
		doTCMET = cms.untracked.bool(False),
		doGenEvent = cms.untracked.bool(False),#put on False when running non-ttbar or when running toptree from reco
		doNPGenEvent = cms.untracked.bool(False),#put on True when running New Physics sample
		doSpinCorrGen = cms.untracked.bool(False),#put on True only if you need SpinCorrelation Variables

		conversionLikelihoodWeightsFile = cms.untracked.string('RecoEgamma/EgammaTools/data/TMVAnalysis_Likelihood.weights.txt'),

		# Draw MC particle tree
		drawMCTree = cms.untracked.bool(False),
		mcTreePrintP4 = cms.untracked.bool(False),
		mcTreePrintPtEtaPhi = cms.untracked.bool(False),
		mcTreePrintVertex = cms.untracked.bool(False),
		mcTreePrintStatus = cms.untracked.bool(False),
		mcTreePrintIndex = cms.untracked.bool(False),
		mcTreeStatus = cms.untracked.vint32( 3 ),	# accepted status codes

	
		# MC particles acceptance cuts
		electronMC_etaMax = cms.double(3.0),
		electronMC_ptMin = cms.double(2.0),
		muonMC_etaMax = cms.double(3.0),
		muonMC_ptMin = cms.double(0.0),
		jetMC_etaMax = cms.double(6.0),
		jetMC_ptMin = cms.double(5.0),
	),

	producersNames = cms.PSet(
		hltProducer1st = cms.InputTag("TriggerResults","","REDIGI"),
		hltProducer2nd = cms.InputTag("TriggerResults","","HLT"),
		hltProducer3rd = cms.InputTag("TriggerResults","","HLTdummy"),
		hltProducer4th = cms.InputTag("TriggerResults","","HLTdummy2"),
    pileUpProducer = cms.InputTag("addPileupInfo"),
		genParticlesProducer = cms.InputTag("genParticles"),
		primaryVertexProducer = cms.InputTag("offlinePrimaryVertices"),
		vcaloJetProducer = cms.untracked.vstring("cleanPatJetsAK5Calo"),
		vgenJetProducer = cms.untracked.vstring("ak5GenJets"),
		vpfJetProducer = cms.untracked.vstring("cleanPatJetsAK5PF"),
		vJPTJetProducer = cms.untracked.vstring("selectedPatJetsAK5JPT"),
		vmuonProducer = cms.untracked.vstring("selectedPatMuons"),
		velectronProducer = cms.untracked.vstring("selectedPatElectrons"),# if electronTriggerMatching == true, change the electron inputTag to "cleanPatElectronsTriggerMatch"
		CalometProducer = cms.InputTag("patMETs"),
    vpfmetProducer = cms.untracked.vstring("patMETsPF"),
    TCmetProducer = cms.InputTag("patMETsTC"),
		genEventProducer = cms.InputTag("genEvt"),
		generalTrackLabel = cms.InputTag("generalTracks"), # to calculate the conversion flag
		electronNewId = cms.untracked.bool(True) #for recent electronID recommanded by EGamma. still Not accepted by Top group
	)
)

process.p = cms.Path(process.analysis)
