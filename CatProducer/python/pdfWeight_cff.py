import FWCore.ParameterSet.Config as cms

pdfWeight = cms.EDProducer("PDFWeightsProducer",
    pdfName = cms.string("cteq66.LHgrid"),
    generatedPdfName = cms.string("cteq66.LHgrid"),
)

