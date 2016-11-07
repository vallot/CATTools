import FWCore.ParameterSet.Config as cms

def customise(process):
    process.out = cms.OutputModule("PoolOutputModule",
        outputCommands = cms.untracked.vstring('drop *',
            'keep *_*_*_CATeX',
            'keep *_partonTop_*_*',
            'keep *_pseudoTop_*_*',
            'keep *_*_genTtbarId_*'),
        fileName = cms.untracked.string("out.root"),
        SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('p'),
        ),
    )

    if hasattr(process, 'outPath'): process.outPath += process.out
    else: process.outPath = cms.EndPath(process.out)
