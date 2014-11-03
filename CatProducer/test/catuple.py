from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("CATTools.CatProducer.catCandidates_cff")

runOnMC=True
postfix = "PFlow"
jetAlgo="AK5"

from Configuration.AlCa.autoCond import autoCond
if runOnMC:
    process.GlobalTag.globaltag = autoCond['startup']
    jecLevels = ['L1FastJet','L2Relative','L3Absolute']
else:
    process.GlobalTag.globaltag = autoCond['com10']
    jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, jetCorrections=("AK5PFchs", jecLevels),
          runOnMC=runOnMC, postfix=postfix, typeIMetCorrections=False)

if not runOnMC:
    removeMCMatchingPF2PAT( process, '' )

process.p = cms.Path(getattr(process,"patPF2PATSequence"+postfix)+process.makeCatCandidates)

from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrectionInPF2PAT
enablePileUpCorrectionInPF2PAT( process, postfix)

from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_cat*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep *_goodOfflinePrimaryVertices*_*_*',
    'keep *_TriggerResults_*_PAT',
#    'keep *_*_*_PAT',
    )

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

#from miniAOD_tools import *
#miniAOD_customizeCommon(process)

process.selectedPatJetsPFlow.cut = cms.string("pt > 20")
process.selectedPatMuonsPFlow.cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))") 
process.selectedPatElectronsPFlow.cut = cms.string("pt > 5") 
process.selectedPatTausPFlow.cut = cms.string("pt > 18. && tauID('decayModeFinding')> 0.5")
process.selectedPatPhotonsPFlow.cut = cms.string("pt > 5")

process.patJetPartonMatchPFlow.mcStatus = [ 3, 23 ]
process.patPFParticlesPFlow.embedGenMatch = cms.bool(True)

process.maxEvents.input = -1
process.out.fileName = 'catTuple.root'

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.catJets.src = cms.InputTag("selectedPatJetsPFlow")
process.catMuons.src = cms.InputTag("selectedPatMuonsPFlow")
process.catElectrons.src = cms.InputTag("selectedPatElectronsPFlow")
process.catPhotons.src = cms.InputTag("patPhotons")
process.catMETs.src = cms.InputTag("patMETsPFlow")


process.source.fileNames = cms.untracked.vstring(
'file:/pnfs/user/kraft_data/FEEEC639-4A98-E211-BE1C-002618943919.root'
#'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/001A0DC8-C313-E211-BCCB-00261894397B.root'
#'file:/pnfs/user/qcd/QCD_Pt-15to3000_Tune1_Flat_8TeV_pythia8_AODSIM_PU_S7_START52_V9-v1/02F7EBD9-A09E-E111-AF77-003048C692C0.root'
    )
#                                         ##
