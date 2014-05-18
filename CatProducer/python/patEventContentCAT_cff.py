import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patTriggerEventContent, patTriggerStandAloneEventContent

patEventContentCAT = cms.untracked.vstring()

patEventContentCAT.extend([
    'keep patJets_selectedPatJetsPF2PAT__PAT',
    'keep patMuons_selectedPatMuonsPF2PAT__PAT',
    'keep patElectrons_selectedPatElectronsPF2PAT__PAT',
    'keep patMETs_patType1CorrectedPFMetPF2PAT__PAT',
    'keep *_selectedPatPhotons*_*_*',
    #'keep *_photons*_*_*',
    'keep double_kt6PFJets_rho_*',
    'keep *_goodOfflinePrimaryVertices*_*_*',
    'keep *_TriggerResults_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep recoGenJets_ak5GenJetsNoNu_*_*',
    'keep *_bFlavorHistoryProducer_*_*',
    'keep *_cFlavorHistoryProducer_*_*',
    'keep *_flavorHistoryFilter_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep recoTracks_generalTracks_*_*',
    'keep recoPFCandidates_selectedPatJetsPF2PAT_pfCandidates_PAT',
    'keep recoGenParticles_genParticles_*_*',
    'keep triggerTriggerEvent_*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep edmMergeableCounter_*_*_*',
    ## we need following output collections below in order to produce top tree photon object from PAT.
    ## it is OK with PAT2TOP_cfg.py without these collections as we can get directly from AOD
    ## in order to avoid saving these, we need to modify PAT photon object
    ## for the time being, let's save it. It should not matter as we don't save PAT.
    'keep *_allConversions_*_*',
    'keep *_offlineBeamSpot_*_*',   
    'keep *_gsfElectron*_*_*', 
    'keep recoPhotonCores_photonCore_*_*',
    'keep *_phPFIsoValueCharged03PFIdPFIso_*_*',
    'keep *_phPFIsoValueGamma03PFIdPFIso_*_*',
    'keep *_phPFIsoValueNeutral03PFIdPFIso_*_*',
    'keep *_phPFIsoValuePU03PFIdPFIso_*_*',
    ])

