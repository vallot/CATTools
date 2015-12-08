import FWCore.ParameterSet.Config as cms

catSecVertexs = cms.EDProducer("CATSecVertexProducer",
    muonSrc = cms.InputTag("slimmedMuons"),
    elecSrc = cms.InputTag("slimmedElectrons"),
    trackSrc = cms.InputTag("generalTracks"),
    secvtxSrc = cms.InputTag("slimmedSecondaryVertices"),
    pfSrc  = cms.InputTag("packedPFCandidates"),
    vertexLabel = cms.InputTag("catVertex"),
    track = cms.PSet(
        minPt = cms.double(1.0),
        maxEta = cms.double(2.5),
        chi2 = cms.double(5.),
        nHit = cms.int32(6),
        signif = cms.double(-5),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        chi2 = cms.double(7.),
        minLxy = cms.double(-100),
        maxLxy = cms.double(100),
        signif = cms.double(-5.0),
    ),
    rawMassMin = cms.double(2),
    rawMassMax = cms.double(4),
    massMin = cms.double(2.80),
    massMax = cms.double(3.40),
    mode = cms.int32(3),
    pfmode = cms.int32(1),  ## pfmode : 0 [ lepton + lepton] ,  1 [ lepton+ Charged Hadron] , 2 [ pdgID cut is not applied. ]
    applyCut = cms.bool(True)
)


