import FWCore.ParameterSet.Config as cms

genWeight = cms.EDProducer("GenWeightsProducer",
    enforceUnitGenWeight = cms.bool(False),
    lheEvent = cms.InputTag("externalLHEProducer"),
    genEventInfo = cms.InputTag("generator"),
    pdfName = cms.string("NNPDF30_nlo_as_0118"),
    generatedPdfName = cms.string("NNPDF30_nlo_as_0118"),
    lheWeightIndex = cms.int32(0),
    genWeightIndex = cms.int32(0),
)

