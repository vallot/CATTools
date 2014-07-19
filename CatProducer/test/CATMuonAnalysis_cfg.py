## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:../prod/CAT.root'
    )
)

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('CATMuonAnalysis.root')
)

process.defaultIso = cms.EDAnalyzer("CATMuonAnalysis",
    # input collection
    src = cms.InputTag("catMuons"),
)

#process.weightIso = cms.EDAnalyzer("CATMuonAnalysis",
#    # input collection
#    src = cms.InputTag("catMuonsWeighted"),
#)


process.p = cms.Path(
            process.defaultIso
#            *process.weightIso
)



