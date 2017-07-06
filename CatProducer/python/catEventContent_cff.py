import FWCore.ParameterSet.Config as cms

catEventContent = cms.untracked.vstring()
catEventContentMC = cms.untracked.vstring()
catEventContentRD = cms.untracked.vstring()
catEventContentTOPMC = cms.untracked.vstring()
catEventContentSecVertexs = cms.untracked.vstring()

catEventContent.extend([
    'drop *',
    'keep *_nEventsTotal_*_*',
    'keep *_nEventsFiltered_*_*',
    'keep *_catMuons_*_*',
    'keep *_catElectrons_*_*',
    # turning off since no one uses photon and taus
    'keep *_catPhotons_*_*',
    #'keep *_catTaus_*_*', 
    'keep *_catJets*_*_*',
    'keep *_catFatJets*_*_*',
    'keep *_catMETs*_*_*',
    'keep *_catVertex_*_*',
    'keep *_catTrigger_*_*',
    'keep edmTriggerResults_TriggerResults__HLT',
    'keep edmTriggerResults_TriggerResults__HLT2',
    'keep patPackedTriggerPrescales_patTrigger__*', ## do we need this collection, where?
    'keep *_lumiMask*_*_*',
    'keep *_fixedGridRhoAll_*_', 'keep *_fixedGridRhoFastjetAll_*_*',
    'keep *_BadChargedCandidateFilter_*_*', 'keep *_BadPFMuonFilter_*_*',
    ])

catEventContentRD.extend([
    'keep *_lumiMask*_*_*',
    ])

catEventContentMC.extend([
    'keep recoGenParticles_prunedGenParticles_*_*',
    'keep *_slimmedGenJets_*_*',
    'keep *_genWeight_*_*',
    'keep *_pileupWeight*_*_*',
    #'keep *_matchGenBHadron_*_*',
    #'keep *_matchGenCHadron_*_*',
    ])

catEventContentTOPMC.extend([
    'keep *_GenTtbarCategories_*_*',
    'keep *_GenTtbarCategories*_*_*',
    'keep *_catGenTops_*_*',
    #'keep *_partonTop_*_*', ## Can be built on the fly from prunedGenParticles
    'keep *_pseudoTop_*_*',
    'keep *_genJetHadronFlavour_*_*',
    ])

catEventContentSecVertexs.extend([
    'keep *_catSecVertexs_*_*',
    ])
