import FWCore.ParameterSet.Config as cms

catEventContent = cms.untracked.vstring()
catEventContentExtended = cms.untracked.vstring()

catEventContent.extend([
    'keep *_catMuons_*_*',
    'keep *_catElectrons_*_*',
    'keep *_catPhotons_*_*',
    'keep *_catJets_*_*',
    'keep *_catMETs_*_*',
    'keep *_catMCParticles_*_*',
    'keep *_catGenTops_*_*',
    ])

catEventContentExtended.extend([
    'keep *_cat*_*_*',
    'keep *_offlineSlimmedPrimaryVertices_*_*',
    'keep *_prunedGenParticles_*_*',
    #'keep patTriggerPaths_patTrigger*_*_*',
    #'keep *_goodOfflinePrimaryVertices*_*_*',
    #'keep GenEventInfoProduct_*_*_*',
    #'keep PileupSummaryInfos_*_*_*',
    #'keep *_selectedPatJets_*_*',
    #'keep *_TriggerResults_*_PAT',
    #'keep *_patTrigger*_*_*',
    #'keep *_*_*_PAT',
    ])

