import FWCore.ParameterSet.Config as cms

ttbarDileptonKin = cms.EDProducer("TTLLKinSolutionProducer",
    solver = cms.PSet(algo = cms.string("Default"),), ## Default dummy algorithm
    leptons = cms.InputTag("eventsTTLL", "leptons"), ## lepton in LorentzVector
    jets = cms.InputTag("eventsTTLL", "jets"), ## jet in LorentzVector
    met = cms.InputTag("eventsTTLL", "met"), ## MET pt in float 
    metphi = cms.InputTag("eventsTTLL", "metphi"), ## MET phi in float
)

