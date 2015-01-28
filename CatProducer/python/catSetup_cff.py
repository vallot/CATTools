import FWCore.ParameterSet.Config as cms

def catSetup(process, runOnMC=True, doSecVertex=True):    
    process.load("CATTools.CatProducer.catCandidates_cff")        

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
        process.load("CATTools.CatProducer.pseudoTop_cfi")
        process.out.outputCommands.append("keep *_pdfWeight_*_*")
        process.out.outputCommands.append("keep *_pileupWeight_*_*")
        process.out.outputCommands.append("keep *_pseudoTop_*_*")
        process.out.outputCommands.append("keep *_partonTop_*_*")

        ## using MEtUncertainties to get lepton shifts
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
        
        ## using objectsUncertaintyTools to get jet shifts
        from PhysicsTools.PatUtils.tools.jmeUncertaintyTools import JetMEtUncertaintyTools, addSmearedJets
        from PhysicsTools.PatUtils.tools.objectsUncertaintyTools import addShiftedJetCollections
        import RecoMET.METProducers.METSigParams_cfi as jetResolutions
        jetSmearFileName='PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'
        jetSmearHistogram='pfJetResolutionMCtoDataCorrLUT'
        varyByNsigmas=1.0
        shiftedParticleSequence = cms.Sequence()
        variations = { "":None, "ResUp":-1., "ResDown":1.  }
        for var in variations.keys():
            jetCollectionToKeep = addSmearedJets(process, catJetsSource,
                                       [ "smeared", catJetsSource, var ],
                                       jetSmearFileName,jetSmearHistogram,jetResolutions,
                                       varyByNsigmas, variations[ var ],
                                       shiftedParticleSequence)
        jetCorrLabelUpToL3='ak4PFL1FastL2L3'
        jetCorrLabelUpToL3Res='ak4PFL1FastL2L3Residual'
        jecUncertaintyFile="PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt"
        jecUncertaintyTag='SubTotalMC'
        shiftedJetsCollections, jetsCollectionsToKeep = addShiftedJetCollections(
            process,catJetsSource,catJetsSource,
            jetCorrLabelUpToL3, jetCorrLabelUpToL3Res,
            jecUncertaintyFile, jecUncertaintyTag,
            varyByNsigmas, shiftedParticleSequence)

            
        process.p += shiftedParticleSequence
        process.catMuons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedMuonsEnDown")
        process.catMuons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedMuonsEnUp")
        process.catElectrons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedElectronsEnDown")
        process.catElectrons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedElectronsEnUp")
        process.catJets.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedJetsEnDown")
        process.catJets.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedJetsEnUp")
        process.catJets.smearedResSrc = cms.InputTag("smearedSlimmedJets")
        process.catJets.smearedResDownSrc = cms.InputTag("smearedSlimmedJetsResDown")
        process.catJets.smearedResUpSrc = cms.InputTag("smearedSlimmedJetsResUp")

        process.out.outputCommands.append("drop *_shifted*_*_*")
        process.out.outputCommands.append("drop *_smeared*_*_*")

    if not runOnMC:
        process.makeCatCandidates.remove(process.catGenJets)
        process.catMuons.runOnMC = cms.bool(False)
        process.catElectrons.runOnMC = cms.bool(False)
        process.catJets.runOnMC = cms.bool(False)

    process.p += process.makeCatCandidates

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        process.makeCatCandidates += process.catSecVertexs

    
