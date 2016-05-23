import FWCore.ParameterSet.Config as cms

catMuons = cms.EDProducer("CATMuonProducer",
    src = cms.InputTag("patMuons"),
    mcLabel = cms.InputTag("genParticles"),
    vertexLabel = cms.InputTag("offlinePrimaryVertices"),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    runOnMC = cms.bool(True)
)

catElectrons = cms.EDProducer("CATElectronProducer",
    src = cms.InputTag("patElectrons"),
    mcLabel = cms.InputTag("genParticles"),
    vertexLabel = cms.InputTag('offlinePrimaryVertices'),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    rhoLabel = cms.InputTag("kt6PFJets", "rho"),
    runOnMC = cms.bool(True),
    electronIDNames = cms.vstring()
)

catJets = cms.EDProducer("CATJetProducer",
    src = cms.InputTag("selectedPatJetsPFlow"),
    shiftedEnDownSrc = cms.InputTag("shiftedSlimmedJetsEnDown"),
    shiftedEnUpSrc = cms.InputTag("shiftedSlimmedJetsEnUp"),
    smearedResSrc = cms.InputTag("smearedSlimmedJets"),
    smearedResDownSrc = cms.InputTag("smearedSlimmedJetsResDown"),
    smearedResUpSrc = cms.InputTag("smearedSlimmedJetsResUp"),
    runOnMC = cms.bool(True),
    btagNames = cms.vstring()
##       #see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance
##       'trackCountingHighPurBJetTags', #0
##       'jetProbabilityBJetTags', #1
##       'combinedSecondaryVertexBJetTags', #2
##  #     'combinedSecondaryVertexV1BJetsTags', #3
##  #     'combinedSecondaryVertexSoftPFLeptonV1BJetTags', #4
##  #     'combinedSecondaryVertexIVFV2BJetTags', #5
## #old for Run 1
##  #     'trackCountingHighEffBJetTags',#0
##  #     'trackCountingHighPurBJetTags',#1
##  #     'jetProbabilityBJetTags',#2
##  #     'jetBProbabilityBJetTags',#3
##  #     'simpleSecondaryVertexHighEffBJetTags',#4
##  #     'simpleSecondaryVertexHighPurBJetTags',#5
##  #     'combinedSecondaryVertexBJetTags',#6
##  #    'combinedSecondaryVertexMVABJetTags'#7
##  #     'softPFMuonBJetTags', #8
##  #     'softPFElectronBJetTags', #9
##  #     'softPFElectronByIP3dBJetTags', #10
##  #     'softPFElectronByPtBJetTags', #11
##  #     'softPFMuonByPtBJetTags', #12
##     ),
)

catPhotons = cms.EDProducer("CATPhotonProducer",
    src = cms.InputTag("selectedPatPhotons"),
)

catTaus = cms.EDProducer("CATTauProducer",
    src = cms.InputTag("selectedPatTaus"),
)

catMETs = cms.EDProducer("CATMETProducer",
    src = cms.InputTag("patMETsPFlow"),
)

catGenJets = cms.EDProducer("CATGenJetProducer",
    src = cms.InputTag("ak5GenJets"),
    pt = cms.double(10),
    eta = cms.double(2.5)
)

catMCParticles = cms.EDProducer("CATMCParticleProducer",
    src = cms.InputTag("genParticles"),
    pt = cms.double(10),
    eta = cms.double(2.5)
)

#catGenTops = cms.EDProducer("CATGenTopProducer",
#    genJetLabel = cms.InputTag("slimmedGenJets"),
#    mcParticleLabel = cms.InputTag("prunedGenParticles"),
#)

catSecVertexs = cms.EDProducer("CATSecVertexProducer",
    muonSrc = cms.InputTag("selectedPatMuonsPFlow"),
    elecSrc = cms.InputTag("selectedPatElectronsPFlow"),
    trackSrc = cms.InputTag("generalTracks"),
    vertexLabel = cms.InputTag("goodOfflinePrimaryVertices"),
    pfmuonSrc = cms.InputTag("pfAllMuonsPFlow"),
    pfelecSrc = cms.InputTag("pfAllElectronsPFlow"),
    track = cms.PSet(
        minPt = cms.double(1.0),
        maxEta = cms.double(2.5),
        ## chi2 = cms.double(5.),
        ## nHit = cms.int32(6),
        ## signif = cms.double(-5),
        ## DCA = cms.double(1.),
        chi2 = cms.double(100),
        nHit = cms.int32(6),
        signif = cms.double(-5),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        ## chi2 = cms.double(7.),
        ## minLxy = cms.double(-100),
        ## maxLxy = cms.double(100),
        ## signif = cms.double(-5.0),
        chi2 = cms.double(100.),
        minLxy = cms.double(-100),
        maxLxy = cms.double(100),
        signif = cms.double(100),
    ),
    rawMassMin = cms.double(1),
    rawMassMax = cms.double(10),
    massMin = cms.double(1),
    massMax = cms.double(10),
)

makeCatCandidates =  cms.Sequence( 
    catMuons*
    catElectrons*
    catJets*
    catMETs
#    catGenJets*
#    catMCParticles
    #catPhotons*
    #catTaus*
    # MC information below
    # dont need for now
    #catGenTops*
    #catSecVertexs
) 
