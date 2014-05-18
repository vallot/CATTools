# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

#this is to test configuration file at CERN
#from CATTools.CatProducer.datasetToSource import *
## This is used to get the correct global tag below, and to find the files
## It is *reset* automatically by ProductionTasks, so you can use it after the ProductionTasksHook
#datasetInfo = ('cmgtools_group', '/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B','.*root')
#process.source = datasetToSource(
#    *datasetInfo
#    )
#process.source.fileNames = process.source.fileNames[:20]

process.source.fileNames = [
#    '/store/mc/Summer12_DR53X/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/0024E066-2BEA-E111-B72F-001BFCDBD11E.root'
#    '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/FED775BD-B8E1-E111-8ED5-003048C69036.root',
    #T2 at Belgium
    '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/4085E811-F197-E211-8A95-002618943953.root',
    #data at CERN
    #'/store/data/Run2012A/DoubleMu/AOD/22Jan2013-v1/30000/FEF469F7-0882-E211-8351-0026189438E6.root'
    #AOD at eos  
#    '/store/caf/user/tjkim/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/70300E2E-27D2-E111-92BD-001E67397AE4.root'
]

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")


runOnMC = True 
runOnFastSim = False 

###############################
####### Global Setup ##########
###############################

process.load("FWCore.Framework.test.cmsExceptionsFatal_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.destinations = ['cerr']
process.MessageLogger.statistics = []
process.MessageLogger.fwkJobReports = []
process.MessageLogger.categories=cms.untracked.vstring('FwkJob'
							,'FwkReport'
							,'FwkSummary'
			                               )

process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options = cms.untracked.PSet(
		                 wantSummary = cms.untracked.bool(True)
	         	 	)

### Set the global tag from the dataset name
if runOnFastSim is False:
  from CATTools.CatProducer.Tools.getGlobalTag import getGlobalTagByDataset
  process.GlobalTag.globaltag = getGlobalTagByDataset( runOnMC, process.source.fileNames[0])
else:
  process.GlobalTag.globaltag = cms.string('START53_V27::All')


##-------------------- Import the Jet RECO modules ----------------------- ## this makes cmsRun crash
##
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

###############################
#### Load MVA electron Id #####
###############################

### UserCode area is not supported anymore ### 
#process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
### Use CMSSW area 
process.load('EgammaAnalysis/ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence( process.mvaTrigV0 + process.mvaNonTrigV0 )

###############################
####### PF2PAT Setup ##########
###############################

# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PF2PAT"
usePF2PAT(process,runPF2PAT=True, jetAlgo="AK5", runOnMC=runOnMC, postfix=postfix, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True)

# TOP projection
process.pfIsolatedMuonsPF2PAT.isolationCut = cms.double(0.2)
process.pfIsolatedMuonsPF2PAT.doDeltaBetaCorrection = True
process.pfSelectedMuonsPF2PAT.cut = cms.string('pt > 10. && abs(eta) < 2.5')
process.pfIsolatedMuonsPF2PAT.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("muPFIsoValueCharged04PF2PAT"))
process.pfIsolatedMuonsPF2PAT.deltaBetaIsolationValueMap = cms.InputTag("muPFIsoValuePU04PF2PAT")
process.pfIsolatedMuonsPF2PAT.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("muPFIsoValueNeutral04PF2PAT"), cms.InputTag("muPFIsoValueGamma04PF2PAT"))
# leptons for top tree: no isolation requirement 
# "pfMuons" is cloned from "pfIsolatedMuons" but the isolation cut is removed in the main PFBRECO sequence
process.patMuonsPF2PAT.pfMuonSource = "pfMuonsPF2PAT"
process.patMuonsPF2PAT.embedCaloMETMuonCorrs = False
process.patMuonsPF2PAT.embedTcMETMuonCorrs= False

print "process.pfIsolatedMuonsPF2PAT.isolationCut -> "+str(process.pfIsolatedMuonsPF2PAT.isolationCut)

# to use GsfElectrons instead of PF electrons
# this will destory the feature of top projection which solves the ambiguity between leptons and jets because
# there will be overlap between non-PF electrons and jets even though top projection is ON!
useGsfElectrons(process,postfix,"03") # to change isolation cone size to 0.3 as it is recommended by EGM POG, use "04" for cone size 0.4

from CATTools.CatProducer.Tools.tools import *
useRecoMuon(process,postfix,"04")

#process.pfIsolatedElectronsPF2PAT.isolationCut = cms.double(0.2)
#process.pfIsolatedElectronsPF2PAT.doDeltaBetaCorrection = False
#process.pfSelectedElectronsPF2PAT.cut = cms.string('pt > 10. && abs(eta) < 2.5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2')
#process.pfIsolatedElectronsPF2PAT.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"))
#process.pfIsolatedElectronsPF2PAT.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT")
#process.pfIsolatedElectronsPF2PAT.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPF2PAT"), cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT"))

#process.patElectronsPF2PAT.isolationValues = cms.PSet(
#    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"),
#    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPF2PAT"),
#    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT"),
#    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPF2PAT"),
#    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT")
#    )

process.patElectronsPF2PAT.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPF2PAT.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePF2PAT.replace( process.patElectronsPF2PAT, process.eidMVASequence * process.patElectronsPF2PAT )

process.patJetCorrFactorsPF2PAT.payload = 'AK5PFchs'
process.patJetCorrFactorsPF2PAT.levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
process.pfPileUpPF2PAT.checkClosestZVertex = False

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = False

#####################################################################################################
#### Clone the PF2PAT sequence for data-driven QCD estimate, and for Stijn's JetMET service work ####
#####################################################################################################

from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
postfixNoLeptonCleaning = 'NoLeptonCleaning'

# just cloning the first sequence, and enabling lepton cleaning 
cloneProcessingSnippet(process, getattr(process, 'patPF2PATSequencePF2PAT'), postfixNoLeptonCleaning)

getattr(process,"pfNoMuonPF2PATNoLeptonCleaning").enable = False
getattr(process,"pfNoElectronPF2PATNoLeptonCleaning").enable = False 
getattr(process,"pfIsolatedMuonsPF2PATNoLeptonCleaning").combinedIsolationCut = cms.double(999999)
getattr(process,"pfIsolatedMuonsPF2PATNoLeptonCleaning").isolationCut = cms.double(999999)
getattr(process,"pfIsolatedElectronsPF2PATNoLeptonCleaning").combinedIsolationCut = cms.double(999999)
getattr(process,"pfIsolatedElectronsPF2PATNoLeptonCleaning").isolationCut = cms.double(999999)

###############################
###### Bare KT 0.6 jets #######
###############################

#from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
# For electron (effective area) isolation
#process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

###############################
### Add AK5GenJetsNoMuNoNu ####
###############################

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoMuNoNu
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

process.ak5GenJetsNoMuNoNu = ak5GenJets.clone( src = cms.InputTag('genParticlesForJetsNoMuNoNu') )
process.ak5GenJetsSeq = cms.Sequence(genParticlesForJetsNoMuNoNu*process.ak5GenJetsNoMuNoNu)

###############################
#### Selections Setup #########
###############################

# AK5 Jets
#   PF
process.selectedPatJetsPF2PAT.cut = cms.string("pt > 10")
process.selectedPatJetsPF2PATNoLeptonCleaning.cut = cms.string("pt > 10")
#process.selectedPatJetsPF2PATNoPFnoPU.cut = cms.string("pt > 10")
#process.selectedPatJetsAK5Calo.cut = cms.string("pt > 15")

# Flavor history stuff
process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")
process.flavorHistoryFilter.pathToSelect = cms.int32(-1)
process.cFlavorHistoryProducer.matchedSrc = cms.InputTag("ak5GenJetsNoNu")
process.bFlavorHistoryProducer.matchedSrc = cms.InputTag("ak5GenJetsNoNu")

process.prePathCounter = cms.EDProducer("EventCountProducer")
process.postPathCounter = cms.EDProducer("EventCountProducer")

process.load('CATTools.CatProducer.eventCleaning.eventCleaning_cff')

### photon sequence ###
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso, setupPFPhotonIso
process.phoIsoSequence = setupPFPhotonIso(process, 'selectedPatPhotons')

process.photonSequence = cms.Sequence (
    process.makePatPhotons+
    process.selectedPatPhotons+
    process.phoIsoSequence
)

# let it run
process.patseq = cms.Sequence(
    process.prePathCounter*
#    process.kt6PFJetsForIsolation*
    process.goodOfflinePrimaryVertices*
#    process.ak5GenJetsSeq*
    process.primaryVertexFilter * #removes events with no good pv (but if cuts to determine good pv change...)
    process.eventCleaning*
    getattr(process,"patPF2PATSequence"+postfix)* # main PF2PAT
#    getattr(process,"patPF2PATSequence"+postfix+postfixNoLeptonCleaning)* # PF2PAT FOR DATA_DRIVEN QCD
#    getattr(process,"patPF2PATSequence"+postfix+postfixNoPFnoPU)* # PF2PAT FOR JETS WITHOUT PFnoPU
#    process.patDefaultSequence*
    process.flavorHistorySeq*
    process.photonSequence
    )


if runOnMC is False:
    process.patseq.remove( process.flavorHistorySeq )
    process.patJetCorrFactorsPF2PAT.levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])

#if runOnFastSim is True:
#normally this option is only for fastsim simulation (and not full sim) but 
# this isn't implemented yet 
if runOnMC is True:
    process.eventCleaning.remove(process.HBHENoiseFilter)


#################
#### ENDPATH ####
#################

nEventsInit = cms.EDProducer("EventCountProducer")

process.p = cms.Path(
    process.patseq+
    process.postPathCounter
    )

process.out.SelectEvents.SelectEvents = cms.vstring('p')

# rename output file
process.out.fileName = "PAT.root"

# process all the events
process.maxEvents.input = 100 #changed

process.options.wantSummary = False
process.out.dropMetaData = cms.untracked.string("DROPPED")

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

from CATTools.CatProducer.patEventContentCAT_cff import patEventContentCAT
process.out.outputCommands = patEventContentCAT  

#uncomment it to save PAT-tuple
process.outpath = cms.EndPath(
#    process.out
    )

###### TOP TREE #######
#where do we use this database information in top tree? (TJ)
#Data measurements from Summer11
#process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
#process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")

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

process.p += process.analysis
