import FWCore.ParameterSet.Config as cms

def catSetup(process, runOnMC=True, doSecVertex=True):    
    process.load("CATTools.CatProducer.catCandidates_cff")        
    process.p += process.makeCatCandidates

    catJetsSource = "slimmedJets"
    catGenJetsSource = "slimmedGenJets"
    catMuonsSource = "slimmedMuons"
    catElectronsSource = "slimmedElectrons"
    catPhotonsSource = "slimmedPhotons"
    catTausSource = "slimmedTaus"
    catMETsSource = "slimmedMETs"
    catVertexSource = "offlineSlimmedPrimaryVertices"
    catMCsource = "prunedGenParticles"
    catBeamSpot = "offlineBeamSpot"
    catRho = "fixedGridRhoAll"

    process.catJets.src = cms.InputTag(catJetsSource)
    process.catMuons.src = cms.InputTag(catMuonsSource)
    process.catMuons.mcLabel = cms.InputTag(catMCsource)
    process.catMuons.vertexLabel = cms.InputTag(catVertexSource)
    process.catMuons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.src = cms.InputTag(catElectronsSource)
    process.catElectrons.vertexLabel = cms.InputTag(catVertexSource)
    process.catElectrons.mcLabel = cms.InputTag(catMCsource)
    process.catElectrons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.rhoLabel = cms.InputTag(catRho)
    process.catPhotons.src = cms.InputTag(catPhotonsSource)
    process.catTaus.src = cms.InputTag(catTausSource)
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catGenJets.src = cms.InputTag(catGenJetsSource)
    process.catSecVertexs.muonSrc = cms.InputTag(catMuonsSource)
    process.catSecVertexs.elecSrc = cms.InputTag(catElectronsSource)
    process.catSecVertexs.vertexLabel = cms.InputTag(catVertexSource)

    process.load("CATTools.CatProducer.recoEventInfo_cfi")
    process.out.outputCommands.append("keep *_recoEventInfo_*_*")
    if runOnMC:
        ## Load MC dependent producers
        process.load("CATTools.CatProducer.pdfWeight_cff")
        process.load("CATTools.CatProducer.pileupWeight_cff")
        process.out.outputCommands.append("keep *_pdfWeight_*_*")
        process.out.outputCommands.append("keep *_pileupWeight_*_*")

        from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
        runType1PFMEtUncertainties(process,
                                   dRjetCleaning = 0.0,
                                    addToPatDefaultSequence=False,
                                    jetCollection=catJetsSource,
                                    electronCollection=catElectronsSource,
                                    photonCollection=catPhotonsSource,
                                    muonCollection=catMuonsSource,
                                    tauCollection=catTausSource,
                                    makeType1p2corrPFMEt=True,
                                    outputModule=None)

        #from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
        #runMEtUncertainties(process,
        #                    electronCollection = cms.InputTag(catElectronsSource),
        #                    jetCollection=cms.InputTag(catJetsSource),
        #                    muonCollection = cms.InputTag(catMuonsSource),
        #                    tauCollection = cms.InputTag(catTausSource) )
        #process.patDefaultSequencePFlow += process.metUncertaintySequence        

    if not runOnMC:
        process.makeCatCandidates.remove(process.catGenJets)
        process.catMuons.runOnMC = cms.bool(False)
        process.catElectrons.runOnMC = cms.bool(False)
        process.catJets.runOnMC = cms.bool(False)

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        process.makeCatCandidates += process.catSecVertexs
