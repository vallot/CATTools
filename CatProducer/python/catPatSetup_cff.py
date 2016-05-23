import FWCore.ParameterSet.Config as cms
## based on patTuple_PF2PAT_cfg

def catPatConfig(process, runOnMC=True, postfix = "PFlow", jetAlgo="AK5", doTriggerSkim=False):
    if runOnMC:
        jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else:
        jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

    ## total event counter
    process.totaEvents = cms.EDProducer("EventCountProducer")
    process.p = cms.Path(process.totaEvents)

    # met cleaning events
    process.load("RecoMET.METFilters.metFilters_cff")
    process.p += process.metFilters
    #process.load("CATTools.CatProducer.eventCleaning.scrapingFilter_cfi")
    #process.p += process.scrapingFilter
        
    # from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorPFnoPU2012
    # change pvCollection to goodOfflinePrimaryVertices
    # and process.pfPileUpPFlow.checkClosestZVertex = False
    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT,removeMCMatchingPF2PAT
    usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, jetCorrections=("AK5PFchs", jecLevels),
              pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
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

    ### skim for qcd data
    if not runOnMC and doTriggerSkim:
        process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
        process.hltHighLevel.HLTPaths = ['HLT_PFJet80_v*','HLT_PFJet140_v*','HLT_PFJet320_v*'] # qcd only
        process.hltHighLevel.andOr = cms.bool(True)
        #process.hltJet.throw = cms.bool(True)
        process.p += process.hltHighLevel
    
    process.p += getattr(process,"patPF2PATSequence"+postfix)
    # temp fix for photons since they are not done with PF2PAT
    process.p += process.photonMatch + process.patPhotons + process.selectedPatPhotons

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
    process.selectedPatTausPFlow.cut = cms.string("pt > 18. && tauID('decayModeFinding')> 0.5")

    process.patJetPartonMatchPFlow.mcStatus = [ 3, 23 ]
    process.patPFParticlesPFlow.embedGenMatch = cms.bool(True)

    ## adding in user variables ## temp to add into objects
    process.patMuonsPFlow.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03PFlow","muPFIsoValueNeutral03PFlow","muPFIsoValueGamma03PFlow","muPFIsoValuePU03PFlow","muPFIsoValueChargedAll03PFlow")

    process.patElectronsPFlow.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFIdPFlow","elPFIsoValueNeutral03PFIdPFlow","elPFIsoValueGamma03PFIdPFlow","elPFIsoValuePU03PFIdPFlow","elPFIsoValueChargedAll03PFIdPFlow")

    process.patJetsPFlow.addTagInfos = cms.bool(True)
    process.patJets.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosAODPFlow"))
    process.patJetsPFlow.userData.userFunctions = cms.vstring(
        "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
        "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0",
        "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
        "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).nTracks() : 0",
        "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
        "tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).value() : 0",
        "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
        "tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).error() : 0",
    )
    process.patJetsPFlow.userData.userFunctionLabels = cms.vstring('vtxMass','vtxNtracks','vtx3DVal','vtx3DSig')

    ## adding pileup jet id
    process.load("CMGTools.External.pujetidsequence_cff")
    process.puJetMva.jets = cms.InputTag("selectedPatJetsPFlow")
    process.puJetId.jets = cms.InputTag("selectedPatJetsPFlow")
    process.puJetIdChs.jets = cms.InputTag("selectedPatJetsPFlow")
    process.puJetIdChs.vertexes = 'goodOfflinePrimaryVertices'
    process.puJetMvaChs.jets = cms.InputTag("selectedPatJetsPFlow")
    process.puJetMvaChs.vertexes = 'goodOfflinePrimaryVertices'
    process.p += process.puJetIdSqeuence + process.puJetIdSqeuenceChs
