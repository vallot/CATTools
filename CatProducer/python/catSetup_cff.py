import FWCore.ParameterSet.Config as cms

def catSetup(process, runOnMC=True, useMiniAOD = True, doSecVertex=True):    
    catJetsSource = "selectedPatJetsPFlow"
    catGenJetsSource = "ak5GenJets"
    catMuonsSource = "selectedPatMuonsPFlow"
    catElectronsSource = "selectedPatElectronsPFlow"
    catPhotonsSource = "selectedPatPhotons"
    catMETsSource = "patMETsPFlow"
    catVertexSource = "offlinePrimaryVertices"
    catMCsource = "genParticles"
    catBeamSpot = "offlineBeamSpot"
    catRho = "kt6PFJets"

    if not useMiniAOD:
        process.load("CATTools.CatProducer.eventCleaning.eventCleaning_cff")
        process.p += process.eventCleaning

    process.load("CATTools.CatProducer.catCandidates_cff")        
    process.p += process.makeCatCandidates

    if useMiniAOD:
        catJetsSource = "slimmedJets"
        catGenJetsSource = "slimmedGenJets"
        catMuonsSource = "slimmedMuons"
        catElectronsSource = "slimmedElectrons"
        catPhotonsSource = "slimmedPhotons"
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
    process.catPhotons.src = cms.InputTag(catPhotonsSource)
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catGenJets.src = cms.InputTag(catGenJetsSource)
    process.catMCParticles.src = cms.InputTag(catMCsource)
    process.catElectrons.rhoLabel = cms.InputTag(catRho)

    if not runOnMC:
        process.makeCatCandidates.remove(process.catGenJets)
        process.makeCatCandidates.remove(process.catMCParticles)
        process.catMuons.runOnMC = cms.bool(False)
        process.catElectrons.runOnMC = cms.bool(False)

    if not doSecVertex:
        process.makeCatCandidates.remove(process.catSecVertexs)

    ## cuts on selected Pat objects
    if not useMiniAOD:
        getattr(process,catJetsSource).cut = cms.string("pt > 20")
        getattr(process,catMuonsSource).cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))") 
        getattr(process,catElectronsSource).cut = cms.string("pt > 5") 
        getattr(process,catPhotonsSource).cut = cms.string("pt > 5")

