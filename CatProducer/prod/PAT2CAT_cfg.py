## import skeleton process
from CATTools.CatProducer.catTemplate_cfg import *
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("CATTools.CatProducer.catCandidates_cff")

from CATTools.CatProducer.Tools.tools import *

useMiniAOD = True
if useMiniAOD:
  miniAOD(process)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
"file:/pnfs/user/jlee/DYJetsToLL_M-50_13TeV-madgraph-pythia8/0432E62A-7A6C-E411-87BB-002590DB92A8.root")
)

#muons with weighted isolation method only for AOD sample
if not useMiniAOD:
  addMuonWeighted(process)

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


from CATTools.CatProducer.catEventContent_cff import *
if useMiniAOD:
  process.out.outputCommands = catEventContent
else:
  process.out.outputCommands = catEventContentExtended
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
