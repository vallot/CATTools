import FWCore.ParameterSet.Config as cms
## based on patTuple_PF2PAT_cfg

def catPatConfig(process, runOnMC=True, postfix = "PFlow", jetAlgo="AK5"):
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC:
        process.GlobalTag.globaltag = autoCond['startup']
        jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else:
        process.GlobalTag.globaltag = autoCond['com10']
        jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

    ## total event counter
    process.totaEvents = cms.EDProducer("EventCountProducer")

    process.p = cms.Path(process.totaEvents)
    
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
    process.out.outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_cat*_*_*',
        'keep *_goodOfflinePrimaryVertices*_*_*',
        'keep GenEventInfoProduct_*_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep edmMergeableCounter_*_*_*',
        'keep patTriggerPaths_patTrigger*_*_*',
        #'keep *_selectedPatJets_*_*',
        #'keep *_TriggerResults_*_PAT',
        #'keep *_patTrigger*_*_*',
        #'keep *_*_*_PAT',
        )
    process.out.fileName = cms.untracked.string('catTuple.root')
