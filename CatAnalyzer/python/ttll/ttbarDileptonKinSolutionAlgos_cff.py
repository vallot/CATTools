import FWCore.ParameterSet.Config as cms

ttbarDileptonKinAlgoPSetCMSKin = cms.PSet(
    algo = cms.string("CMSKin"),
    tMassBegin = cms.double(100),
    tMassEnd   = cms.double(300),
    tMassStep  = cms.double(0.5),
    #nuPars = cms.vdouble(18.41,61.61,23.68,61.64,23.65), # parton 13TeV 25ns powheg, produced by Youn
    nuPars = cms.vdouble(27.23,53.88,19.92,53.89,19.9), # particle 13TeV 50ns powheg, produced by Youn
    #nuPars = cms.vdouble(30.641,57.941,22.344,57.533,22.232), # 7?8?TeV from DESY
    #nuPars = cms.vdouble(30.7137,56.2880,23.0744,59.1015,24.9145), # 7TeV pythia from CMSSW default
)

ttbarDileptonKinAlgoPSetMT2 = cms.PSet(
    algo = cms.string("MT2"),
)

ttbarDileptonKinAlgoPSetMAOS = cms.PSet(
    solver = cms.string("MAOS"),
)

ttbarDileptonKinAlgoPSetDESYSmeared = cms.PSet(
    algo = cms.string("DESYSmeared"),
    #inputTemplatePath = cms.string("CATTools/CatAnalyzer/data/desyKinRecoInput.root"),
    inputTemplatePath = cms.string("CATTools/CatAnalyzer/data/KoreaDesyKinRecoInput.root"),
    nTrial = cms.int32(100),
    maxLBMass = cms.double(360),
    mTopInput = cms.double(172.5),
)

ttbarDileptonKinAlgoPSetDESYMassLoop = cms.PSet(
    algo = cms.string("DESYMassLoop"),
    tMassBegin = cms.double(100),
    tMassEnd   = cms.double(300),
    tMassStep  = cms.double(0.5),
)

