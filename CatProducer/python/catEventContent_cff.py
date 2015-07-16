import FWCore.ParameterSet.Config as cms

catEventContent = cms.untracked.vstring()
catEventContentMC = cms.untracked.vstring()
catEventContentAODMC = cms.untracked.vstring()
catEventContentSecVertexs = cms.untracked.vstring()

catEventContent.extend([
    'drop *',
    'keep *_catMuons_*_*',
    'keep *_catElectrons_*_*',
    'keep *_catPhotons_*_*',
    'keep *_catJets_*_*',
    'keep *_catMETs_*_*',
    'keep recoVertexs_offlineSlimmedPrimaryVertices_*_*',
    'keep *_recoEventInfo_*_*',
    'drop *_shifted*_*_*',
    'drop *_smeared*_*_*',
    ])

catEventContentMC.extend([
    'keep recoGenParticles_prunedGenParticles_*_*',
    'keep *_slimmedGenJets_*_*',
    'keep *_pdfWeight_*_*',
    'keep *_pileupWeight_*_*',
    #'keep *_pseudoTop_*_*',
    #'keep *_partonTop_*_*',
    ])

catEventContentAODMC.extend([
    'keep *_catGenTops_*_*',
    ])

catEventContentSecVertexs.extend([
    'keep *_catSecVertexs_*_*',
    ])
