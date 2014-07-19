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
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_0_6_patch1/RelValZMM_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/862FBA61-8702-E411-8EE8-003048D25BA6.root')
)

## for muon isolation study with weighted method
print "4508 PR needs to be merged"
print "git cms-merge-topic 4508"
 
process.load("CommonTools.ParticleFlow.deltaBetaWeights_cff")

from PhysicsTools.PatAlgos.tools.helpers import loadWithPostfix
loadWithPostfix(process,'RecoMuon.MuonIsolation.muonPFIsolation_cff',"Weighted")

process.patMuonsWeighted = process.patMuons.clone()
process.catMuonsWeighted = process.catMuons.clone()
process.catMuonsWeighted.src = 'patMuonsWeighted'

process.muPFIsoDepositNeutralWeighted.ExtractorPSet.inputCandView = 'pfWeightedNeutralHadrons'
process.muPFIsoDepositGammaWeighted.ExtractorPSet.inputCandView = 'pfWeightedPhotons'

process.patMuonsWeighted.isoDeposits = cms.PSet(
    pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedWeighted" ),
    pfChargedAll = cms.InputTag("muPFIsoDepositChargedAllWeighted" ),
    pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPUWeighted" ),
    pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralWeighted" ),
    pfPhotons = cms.InputTag("muPFIsoDepositGammaWeighted" ),
    )

process.patMuonsWeighted.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04Weighted"),
    pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04Weighted"),
    pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04Weighted" ),
    pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04Weighted" ),
    pfPhotons = cms.InputTag("muPFIsoValueGamma04Weighted" ),
    )

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

process.out.outputCommands = ['keep *_cat*_*_*']
