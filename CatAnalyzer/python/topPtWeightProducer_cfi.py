import FWCore.ParameterSet.Config as cms

topPtWeight = cms.EDProducer("TopPtWeightProducer",
    src = cms.InputTag("partonTop"),
    ## parameters for the pt weighting: sqrt(exp(params[0]+params[1]*pt1)*exp(params[0]+params[1]*pt2))
    ## NOTE: only 2 parameters are allowed at this moment
    paramsAny = cms.vdouble(0.156, -0.00137),
    paramsLL = cms.vdouble(0.148, -0.00129),
    paramsLJ = cms.vdouble(0.159, -0.00141),
)
