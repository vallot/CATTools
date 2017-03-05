import FWCore.ParameterSet.Config as cms

topPtWeight = cms.EDProducer("TopPtWeightProducer",
    src = cms.InputTag("partonTop"),
    ## parameters for the pt weighting: sqrt(exp(params[0]+params[1]*pt1)*exp(params[0]+params[1]*pt2))
    ## NOTE: only 2 parameters are allowed at this moment
    ## Parameters for 13TeV
    paramsAny = cms.vdouble(0.0615, -0.0005),
    paramsLJ  = cms.vdouble(0.0615, -0.0005),
    paramsLL  = cms.vdouble(0.0615, -0.0005),
    ## Parameters for 8TeV
    #paramsAny = cms.vdouble(0.156, -0.00137),
    #paramsLJ  = cms.vdouble(0.159, -0.00141),
    #paramsLL  = cms.vdouble(0.148, -0.00129),
    ## Parameters for 7TeV
    #paramsAny = cms.vdouble(0.199, -0.00166),
    #paramsLJ  = cms.vdouble(0.174, -0.00137),
    #paramsLL  = cms.vdouble(0.222, -0.00197),
)
