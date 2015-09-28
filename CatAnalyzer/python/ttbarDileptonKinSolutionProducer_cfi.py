import FWCore.ParameterSet.Config as cms

ttbarDileptonKin = cms.EDProducer("TTbarDileptonKinSolutionProducer",
    solver = cms.PSet(algo = cms.string("Default"),), ## Default dummy algorithm
    leptons = cms.InputTag("leptons"), ## lepton in CandidatePtr
    jets = cms.InputTag("jets"), ## jet in CandidatePtr
    met = cms.InputTag("met"), ## MET pt in float 
    metphi = cms.InputTag("metphi"), ## MET phi in float
)

