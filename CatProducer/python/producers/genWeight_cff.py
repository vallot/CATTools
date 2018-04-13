import FWCore.ParameterSet.Config as cms

genWeight = cms.EDProducer("GenWeightsProducer",
    doLOPDFReweight = cms.bool(False),
    enforceUnitGenWeight = cms.bool(False),
    keepFirstOnly = cms.bool(True),
    lheEvent = cms.InputTag("externalLHEProducer"),
    genEventInfo = cms.InputTag("generator"),
    pdfName = cms.string("NNPDF31_nnlo_as_0118"),
    generatedPdfName = cms.string("NNPDF31_nnlo_as_0118"),
)

