import FWCore.ParameterSet.Config as cms

def catSetup(process, runOnMC=True, doSecVertex=True):    
    process.load("CATTools.CatProducer.eventCleaning.eventCleaning_cff")
    process.load("CATTools.CatProducer.catCandidates_cff")
    process.p += process.eventCleaning + process.makeCatCandidates

    catJetsSource = "selectedPatJetsPFlow"
    catGenJetsSource = "ak5GenJets"
    catMuonsSource = "selectedPatMuonsPFlow"
    catElectronsSource = "selectedPatElectronsPFlow"
    catPhotonsSource = "selectedPatPhotons"
    catTausSource = "selectedPatTauPFlow"
    catMETsSource = "patMETsPFlow"
    catVertexSource = "offlinePrimaryVertices"
    catMCsource = "genParticles"
    catBeamSpot = "offlineBeamSpot"
    #catRho = "kt6PFJets"

    process.catJets.src = cms.InputTag(catJetsSource)
    process.catMuons.src = cms.InputTag(catMuonsSource)
    process.catMuons.mcLabel = cms.InputTag(catMCsource)
    process.catMuons.vertexLabel = cms.InputTag(catVertexSource)
    process.catMuons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.src = cms.InputTag(catElectronsSource)
    process.catElectrons.vertexLabel = cms.InputTag(catVertexSource)
    process.catElectrons.mcLabel = cms.InputTag(catMCsource)
    process.catElectrons.beamLineSrc = cms.InputTag(catBeamSpot)
    #process.catElectrons.rhoLabel = cms.InputTag(catRho)
    process.catPhotons.src = cms.InputTag(catPhotonsSource)
    process.catTaus.src = cms.InputTag(catTausSource)
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catGenJets.src = cms.InputTag(catGenJetsSource)
    process.catSecVertexs.muonSrc = cms.InputTag("pfMuonsFromVertexPFlow")
    process.catSecVertexs.elecSrc = cms.InputTag(catElectronsSource)
    process.catSecVertexs.vertexLabel = cms.InputTag(catVertexSource)

    if not runOnMC:
        process.makeCatCandidates.remove(process.catGenJets)
        process.makeCatCandidates.remove(process.catMCParticles)
        process.catMuons.runOnMC = cms.bool(False)
        process.catElectrons.runOnMC = cms.bool(False)

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        process.makeCatCandidates += process.catSecVertexs
        
    ## cuts on selected Pat objects
    getattr(process,catJetsSource).cut = cms.string("pt > 20")
    process.pfSelectedMuonsPFlow.cut = cms.string("")

    #getattr(process,catMuonsSource).cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))") 
    #getattr(process,catElectronsSource).cut = cms.string("pt > 5") 
    #getattr(process,catPhotonsSource).cut = cms.string("pt > 5")
