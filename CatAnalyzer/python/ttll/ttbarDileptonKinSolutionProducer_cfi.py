import FWCore.ParameterSet.Config as cms

ttbarDileptonKin = cms.EDProducer("TTLLKinSolutionProducer",
    solver = cms.PSet(algo = cms.string("Default"),), ## Default dummy algorithm
    leptons = cms.InputTag("leptons"), ## lepton in LorentzVector
    jets = cms.InputTag("jets"), ## jet in LorentzVector
    met = cms.InputTag("met"), ## MET pt in float 
    metphi = cms.InputTag("metphi"), ## MET phi in float
)

