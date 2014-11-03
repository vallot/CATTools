import FWCore.ParameterSet.Config as cms

catJets = cms.EDProducer("CATJetProducer",
    # input collection
    src = cms.InputTag("patJets"),
    btagType = cms.vstring(
      #see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance
      'trackCountingHighPurBJetTags', #0
      'jetProbabilityBJetTags', #1
      'combinedSecondaryVertexBJetTags', #2
 #     'combinedSecondaryVertexV1BJetsTags', #3
 #     'combinedSecondaryVertexSoftPFLeptonV1BJetTags', #4
 #     'combinedSecondaryVertexIVFV2BJetTags', #5
#old for Run 1
 #     'trackCountingHighEffBJetTags',#0
 #     'trackCountingHighPurBJetTags',#1
 #     'jetProbabilityBJetTags',#2
 #     'jetBProbabilityBJetTags',#3
 #     'simpleSecondaryVertexHighEffBJetTags',#4
 #     'simpleSecondaryVertexHighPurBJetTags',#5
 #     'combinedSecondaryVertexBJetTags',#6
 #    'combinedSecondaryVertexMVABJetTags'#7
 #     'softPFMuonBJetTags', #8
 #     'softPFElectronBJetTags', #9
 #     'softPFElectronByIP3dBJetTags', #10
 #     'softPFElectronByPtBJetTags', #11
 #     'softPFMuonByPtBJetTags', #12
    ),
)

