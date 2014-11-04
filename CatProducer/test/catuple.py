from PhysicsTools.PatAlgos.patTemplate_cfg import *

runOnMC=True
postfix = "PFlow"
jetAlgo="AK5"

process.load("CATTools.CatProducer.eventCleaning.eventCleaning_cff")
process.load("CATTools.CatProducer.catCandidates_cff")
process.totaEvents   = cms.EDProducer("EventCountProducer")

from Configuration.AlCa.autoCond import autoCond
if runOnMC:
    process.GlobalTag.globaltag = autoCond['startup']
    jecLevels = ['L1FastJet','L2Relative','L3Absolute']
else:
    process.GlobalTag.globaltag = autoCond['com10']
    jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, jetCorrections=("AK5PFchs", jecLevels),
          runOnMC=runOnMC, postfix=postfix, typeIMetCorrections=True)

if not runOnMC:
    removeMCMatchingPF2PAT( process, '' )

from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrectionInPF2PAT
enablePileUpCorrectionInPF2PAT( process, postfix, sequence = "patPF2PATSequence"+postfix)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, sequence = "patPF2PATSequence"+postfix )

process.p = cms.Path(process.totaEvents
                     +getattr(process,"patPF2PATSequence"+postfix)
                     +process.photonMatch+process.patPhotons+process.selectedPatPhotons
                     +process.eventCleaning+process.makeCatCandidates)

process.out.outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_cat*_*_*',
    'keep *_goodOfflinePrimaryVertices*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep edmMergeableCounter_*_*_*',
    'keep patTriggerPaths_patTrigger*_*_*',
    #'keep *_selectedPatJets_*_*',
    #'keep *_TriggerResults_*_PAT',
    #'keep *_patTrigger*_*_*',
    #'keep *_*_*_PAT',
    )
process.out.fileName = cms.untracked.string('catTuple.root')

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True
# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False
# enable delta beta correction for muon selection in PF2PAT?
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False

process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
process.patDefaultSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )

process.selectedPatJetsPFlow.cut = cms.string("pt > 20")
process.selectedPatMuonsPFlow.cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))") 
process.selectedPatElectronsPFlow.cut = cms.string("pt > 5") 
process.selectedPatTausPFlow.cut = cms.string("pt > 18. && tauID('decayModeFinding')> 0.5")
process.selectedPatPhotonsPFlow.cut = cms.string("pt > 5")

process.patJetPartonMatchPFlow.mcStatus = [ 3, 23 ]
process.patPFParticlesPFlow.embedGenMatch = cms.bool(True)

process.patMuonsPFlow.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03PFlow","muPFIsoValueNeutral03PFlow","muPFIsoValueGamma03PFlow","muPFIsoValuePU03PFlow","muPFIsoValueChargedAll03PFlow")
process.patElectronsPFlow.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFIdPFlow","elPFIsoValueNeutral03PFIdPFlow","elPFIsoValueGamma03PFIdPFlow","elPFIsoValuePU03PFIdPFlow","elPFIsoValueChargedAll03PFIdPFlow")

process.patJetsPFlow.addTagInfos = cms.bool(True)
process.patJetsPFlow.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0",
"? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).value() : 0",
"? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).error() : 0",
)
process.patJetsPFlow.userData.userFunctionLabels = cms.vstring('secvtxMass','Lxy','LxyErr')

process.catJets.src = cms.InputTag("selectedPatJetsPFlow")
process.catMuons.src = cms.InputTag("selectedPatMuonsPFlow")
process.catElectrons.src = cms.InputTag("selectedPatElectronsPFlow")
process.catPhotons.src = cms.InputTag("selectedPatPhotons")
process.catMETs.src = cms.InputTag("patMETsPFlow")

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.maxEvents.input = 100
process.source.fileNames = cms.untracked.vstring(
'file:/pnfs/user/kraft_data/FEEEC639-4A98-E211-BE1C-002618943919.root',
#'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/001A0DC8-C313-E211-BCCB-00261894397B.root'
#'file:/pnfs/user/qcd/QCD_Pt-15to3000_Tune1_Flat_8TeV_pythia8_AODSIM_PU_S7_START52_V9-v1/02F7EBD9-A09E-E111-AF77-003048C692C0.root'
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/003F2216-14E1-E111-960B-003048C693EE.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/020431EA-28E1-E111-B959-0030487E4ED5.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/02B56FF3-11E1-E111-ABE0-002590494C94.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/02F44798-2DE1-E111-83F3-00266CF327C4.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/088B3B16-38E1-E111-99E6-0030487D5DB1.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/08BEF948-28E1-E111-AAF6-0025901D490C.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/0A2B5181-33E1-E111-8A7C-002590494C74.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/0C3D2D22-43E1-E111-A107-00266CF9C1AC.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/0C608324-36E1-E111-8708-003048F0E5A4.root',
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_Tune4C_Flat_8TeV_pythia8/14085CA3-3DE1-E111-BB95-00266CF270A8.root'
)
