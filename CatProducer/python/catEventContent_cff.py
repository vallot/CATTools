import FWCore.ParameterSet.Config as cms

catEventContent = cms.untracked.vstring()
catEventContentMC = cms.untracked.vstring()
catEventContentAODMC = cms.untracked.vstring()
catEventContentSecVertexs = cms.untracked.vstring()

catEventContent.extend([
    'drop *',
    'keep *_nEventsTotal_*_*',
    'keep *_catMuons_*_*',
    'keep *_catElectrons_*_*',
    'keep *_catPhotons_*_*',
    'keep *_catJets*_*_*',
    'keep *_catMETs*_*_*',
    'keep *_catVertex_*_*',
    #'keep recoVertexs_offlineSlimmedPrimaryVertices_*_*',
    'keep *_catTrigger_*_*',
    'drop *_shifted*_*_*',
    'drop *_smeared*_*_*',
    ])

catEventContentMC.extend([
    'keep recoGenParticles_prunedGenParticles_*_*',
    'keep *_slimmedGenJets_*_*',
    'keep *_pdfWeight_*_*',
    'keep *_pileupWeight_*_*',
    #'keep *_matchGenBHadron_*_*',
    #'keep *_matchGenCHadron_*_*',
    'keep *_GenTtbarCategories_*_*',
    'keep *_pseudoTop_*_*',
    #'keep *_partonTop_*_*',
    ])

catEventContentAODMC.extend([
    'keep *_catGenTops_*_*',
    ])

catEventContentSecVertexs.extend([
    'keep *_catSecVertexs_*_*',
    ])
