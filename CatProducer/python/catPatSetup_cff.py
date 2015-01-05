import FWCore.ParameterSet.Config as cms
## based on patTuple_PF2PAT_cfg

def catPatConfig(process, runOnMC=True, postfix = "PFlow", jetAlgo="AK5"):
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC:
        process.GlobalTag.globaltag = autoCond['startup']
        jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else:
        process.GlobalTag.globaltag = autoCond['com10']
        jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT,removeMCMatchingPF2PAT
    usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, jetCorrections=("AK5PFchs", jecLevels),
            runOnMC=runOnMC, postfix=postfix, typeIMetCorrections=True)

    ## pile up corrections
    from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrectionInPF2PAT
    enablePileUpCorrectionInPF2PAT( process, postfix, sequence = "patPF2PATSequence"+postfix)

    ## electron ID tool
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )

    ## adding trigger info
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
    switchOnTrigger( process, sequence = "patPF2PATSequence"+postfix )

    ## total event counter
    process.totaEvents = cms.EDProducer("EventCountProducer")

    process.p = cms.Path(process.totaEvents
        + getattr(process,"patPF2PATSequence"+postfix)
        # temp fix for photons since they are not done with PF2PAT
        + process.photonMatch + process.patPhotons + process.selectedPatPhotons
    )
    if not runOnMC:
        removeMCMatchingPF2PAT( process, postfix=postfix )
        process.p.remove(process.photonMatch)
    
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
    

    # no taus for now...
    #process.selectedPatTausPFlow.cut = cms.string("pt > 18. && tauID('decayModeFinding')> 0.5")

    process.patJetPartonMatchPFlow.mcStatus = [ 3, 23 ]
    process.patPFParticlesPFlow.embedGenMatch = cms.bool(True)

    ## adding in user variables ## temp to add into objects
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


    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
    process.out.outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_cat*_*_*',
        'keep *_goodOfflinePrimaryVertices*_*_*',
        'keep GenEventInfoProduct_*_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep edmMergeableCounter_*_*_*',
        'keep patTriggerPaths_patTrigger*_*_*',
        'keep recoGenParticles_genParticles__SIM',
        #'keep *_selectedPatJets_*_*',
        #'keep *_TriggerResults_*_PAT',
        #'keep *_patTrigger*_*_*',
        #'keep *_*_*_PAT',
        )
    process.out.fileName = cms.untracked.string('catTuple.root')
