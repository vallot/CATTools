import FWCore.ParameterSet.Config as cms

genWeight = cms.EDProducer("GenWeightsProducer",
    doLOPDFReweight = cms.bool(False),
    enforceUnitGenWeight = cms.bool(False),
    lheEvent = cms.InputTag("externalLHEProducer"),
    genEventInfo = cms.InputTag("generator"),
    pdfName = cms.string("NNPDF30_nlo_as_0118"),
    generatedPdfName = cms.string("NNPDF30_nlo_as_0118"),
)

