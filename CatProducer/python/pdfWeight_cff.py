import FWCore.ParameterSet.Config as cms

pdfWeight = cms.EDProducer("PDFWeightsProducer",
    genEventInfo = cms.InputTag("generator"),
    pdfName = cms.string("NNPDF30_nlo_as_0118"),
    generatedPdfName = cms.string("NNPDF30_nlo_as_0118"),
)

