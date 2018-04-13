import FWCore.ParameterSet.Config as cms

deepcsvWeights = cms.EDProducer("CSVWeightProducer",
    src = cms.InputTag("catJets"),
    leptons = cms.InputTag(""), ## Selected leptons to apply overlap removal
    minPt = cms.double(20),
    maxEta = cms.double(2.41),
    csvSFHF = cms.string("Deepcsv_rwt_fit_hf_v2_final_2018_2_12test.root"),
    csvSFLF = cms.string("Deepcsv_rwt_fit_lf_v2_final_2018_2_12test.root"),
)
