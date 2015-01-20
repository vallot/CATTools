import FWCore.ParameterSet.Config as cms

pdfWeight = cms.EDProducer("PDFWeightsProducer",
    genEventInfo = cms.InputTag("generator"),
    pdfName = cms.string("cteq66.LHgrid"),
    generatedPdfName = cms.string("cteq66.LHgrid"),
)

