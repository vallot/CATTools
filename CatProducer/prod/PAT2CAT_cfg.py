## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("CATTools.CatProducer.catCandidates_cff")

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
#CERN
#       '/store/relval/CMSSW_7_0_6_patch3/RelValZMM_13/GEN-SIM-RECO/PUpmx50ns_PLS170_V6AN1-v2/00000/1EF0EB3F-B412-E411-A7EC-0025905A612A.root'
#        '/store/relval/CMSSW_7_0_7/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PLS170_V7AN1-v1/00000/46E6309D-9516-E411-A4FC-0025905A48F0.root'
    #Kisti
#         'file:/cms/home/tjkim/store/relval/CMSSW_7_0_7/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PLS170_V7AN1-v1/00000/56517256-CB15-E411-9C56-0025905A48F2.root'
'file:/cms/home/jlee/scratch/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp/001A0DC8-C313-E211-BCCB-00261894397B.root'
    )
)

### add different cone size
#let's use userIsolation function to use different cone size for the time being
#userIsolation("pat::User1Iso") for chargedHadronIso()
#userIsolation("pat::User2Iso") for neutralHadronIso()
#userIsolation("pat::User3Iso") for photonIso()
#userIsolation("pat::User4Iso") for puChargedHadronIso()
#userIsolation("pat::User5Iso") for particleIso()
#for muon
## process.patMuons.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03","muPFIsoValueNeutral03","muPFIsoValueGamma03","muPFIsoValuePU03","muPFIsoValueChargedAll03")
## #for electron
## process.patElectrons.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFId","elPFIsoValueNeutral03PFId","elPFIsoValueGamma03PFId","elPFIsoValuePU03PFId","elPFIsoValueChargedAll03PFId")
## process.patElectrons.isolationValuesNoPFId.user = cms.VInputTag("elPFIsoValueCharged03NoPFId","elPFIsoValueNeutral03NoPFId","elPFIsoValueGamma03NoPFId","elPFIsoValuePU03NoPFId","elPFIsoValueChargedAll03NoPFId")
## ###

## process.load("CommonTools.ParticleFlow.deltaBetaWeights_cff")

## from PhysicsTools.PatAlgos.tools.helpers import loadWithPostfix
## loadWithPostfix(process,'RecoMuon.MuonIsolation.muonPFIsolation_cff',"Weighted")

## process.patMuonsWeighted = process.patMuons.clone()
## process.catMuonsWeighted = process.catMuons.clone()
## process.catMuonsWeighted.src = 'patMuonsWeighted'

## process.muPFIsoDepositNeutralWeighted.ExtractorPSet.inputCandView = 'pfWeightedNeutralHadrons'
## process.muPFIsoDepositGammaWeighted.ExtractorPSet.inputCandView = 'pfWeightedPhotons'

## process.patMuonsWeighted.isoDeposits = cms.PSet(
##     pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedWeighted" ),
##     pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralWeighted" ),
##     pfPhotons = cms.InputTag("muPFIsoDepositGammaWeighted" ),
##     pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPUWeighted" ),
##     pfChargedAll = cms.InputTag("muPFIsoDepositChargedAllWeighted" ),
##     )

## process.patMuonsWeighted.isolationValues = cms.PSet(
##     pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04Weighted"),
##     pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04Weighted" ),
##     pfPhotons = cms.InputTag("muPFIsoValueGamma04Weighted" ),
##     pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04Weighted" ),
##     pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04Weighted"),
##     user = cms.VInputTag("muPFIsoValueCharged03Weighted","muPFIsoValueNeutral03Weighted","muPFIsoValueGamma03Weighted","muPFIsoValuePU03Weighted","muPFIsoValueChargedAll03Weighted"),
##     )



#we need following lines for the time being for new btags 
#process.load('CondCore.DBCommon.CondDBSetup_cfi')
#process.BTauMVAJetTagComputerRecord = cms.ESSource('PoolDBESSource',
#    process.CondDBSetup,
#    timetype = cms.string('runnumber'),
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('BTauGenericMVAJetTagComputerRcd'),
#        tag = cms.string('MVAComputerContainer_Retrained53X_JetTags_v2')
#    )),
#    connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
#    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
#)
#process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer('PoolDBESSource','BTauMVAJetTagComputerRecord')

process.patJets.discriminatorSources = [
    cms.InputTag("trackCountingHighPurBJetTags"),
    cms.InputTag("jetProbabilityBJetTags"),
    cms.InputTag("combinedSecondaryVertexBJetTags"),
#    cms.InputTag("combinedSecondaryVertexV1BJetsTags"),
#    cms.InputTag("combinedSecondaryVertexSoftPFLeptonV1BJetTags"),
#    cms.InputTag("combinedSecondaryVertexIVFV2BJetTags")
]
#add b-tag information
process.patJets.addTagInfos = False #let's set False for the time being
process.patJets.tagInfoSources = cms.VInputTag(
  cms.InputTag("secondaryVertexTagInfos")
)
process.patJets.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0",
"? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).value() : 0",
"? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).error() : 0",
)
process.patJets.userData.userFunctionLabels = cms.vstring('secvtxMass','Lxy','LxyErr')


##
## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
#from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
#process.source.fileNames = filesRelValProdTTbarAODSIM
#                                         ##
process.maxEvents.input = -1
#                                         ##
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#                                         ##
process.out.fileName = 'CAT.root'
#                                         ##
#   process.options.wantSummary = False   ##  (to suppress the long output at the end of the job)

process.out.outputCommands = ['keep *_cat*_*_*',
                              'keep *_goodOfflinePrimaryVertices_*_*'
                             ]
