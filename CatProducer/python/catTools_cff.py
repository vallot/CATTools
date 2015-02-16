import FWCore.ParameterSet.Config as cms

def catTool(process, runOnMC=True, doSecVertex=True, useMiniAOD = True):
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

    process.load("CATTools.CatProducer.catCandidates_cff")        
    process.load("CATTools.CatProducer.recoEventInfo_cfi")

    if runOnMC:## Load MC dependent producers
        process.load("CATTools.CatProducer.pdfWeight_cff")
        process.load("CATTools.CatProducer.pileupWeight_cff")
        process.load("CATTools.CatProducer.pseudoTop_cfi")
        if not useMiniAOD:
            process.load("CATTools.CatProducer.genTopProducer_cfi")

    process.catJets.src = cms.InputTag(catJetsSource)
    process.catJets.genJetMatch = cms.InputTag("patJetGenJetMatch")
    process.catTaus.src = cms.InputTag(catTausSource)
    process.catTaus.genJetMatch = cms.InputTag("tauGenJetMatch")
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
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catSecVertexs.muonSrc = cms.InputTag(catMuonsSource)
    process.catSecVertexs.elecSrc = cms.InputTag(catElectronsSource)
    process.catSecVertexs.vertexLabel = cms.InputTag(catVertexSource)

    if runOnMC:
        ## using MEtUncertainties to get lepton shifts
        ## Need to update - Jet calculation was old, would most likely be the same for leptons
        from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
        runType1PFMEtUncertainties(process,
                                    addToPatDefaultSequence=False,
                                    jetCollection=catJetsSource,
                                    electronCollection=catElectronsSource,
                                    photonCollection=catPhotonsSource,
                                    muonCollection=catMuonsSource,
                                    tauCollection=catTausSource,
                                    makeType1p2corrPFMEt=True,
                                    outputModule=None)
        
        process.catMuons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedMuonsEnDown")
        process.catMuons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedMuonsEnUp")
        process.catElectrons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedElectronsEnDown")
        process.catElectrons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedElectronsEnUp")

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        #process.makeCatCandidates += process.catSecVertexs

    
