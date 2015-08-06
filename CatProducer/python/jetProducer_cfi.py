import FWCore.ParameterSet.Config as cms

catJets = cms.EDProducer("CATJetProducer",
    src = cms.InputTag("slimmedJets"),
    btagNames = cms.vstring(),
    ## payloadName should be AK4PFchs, but PHYS14_25_V2 does not have uncertainty 
    payloadName = cms.string(""), 
##       #see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance
##       'trackCountingHighPurBJetTags', #0
##       'jetProbabilityBJetTags', #1
##       'combinedSecondaryVertexBJetTags', #2
##  #     'combinedSecondaryVertexV1BJetsTags', #3
##  #     'combinedSecondaryVertexSoftPFLeptonV1BJetTags', #4
##  #     'combinedSecondaryVertexIVFV2BJetTags', #5
## #old for Run 1
##  #     'trackCountingHighEffBJetTags',#0
##  #     'trackCountingHighPurBJetTags',#1
##  #     'jetProbabilityBJetTags',#2
##  #     'jetBProbabilityBJetTags',#3
##  #     'simpleSecondaryVertexHighEffBJetTags',#4
##  #     'simpleSecondaryVertexHighPurBJetTags',#5
##  #     'combinedSecondaryVertexBJetTags',#6
##  #    'combinedSecondaryVertexMVABJetTags'#7
##  #     'softPFMuonBJetTags', #8
##  #     'softPFElectronBJetTags', #9
##  #     'softPFElectronByIP3dBJetTags', #10
##  #     'softPFElectronByPtBJetTags', #11
##  #     'softPFMuonByPtBJetTags', #12
##     ),
)

catJetsPuppi = cms.EDProducer("CATJetProducer",
    src = cms.InputTag("slimmedJetsPuppi."),
    btagNames = cms.vstring(),
    payloadName = cms.string(""), 
)

#There is a CMS rule that we are supposed to use one module per cfi file.  

