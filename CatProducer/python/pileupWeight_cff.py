import FWCore.ParameterSet.Config as cms
#from   FWCore.PythonUtilities.LumiList import LumiList
from   CATTools.CatProducer.pileupWeight.pileupWeight2018_cff import pileupWeightMap as pileupWeight2018
pileupWeightMap = {}
pileupWeightMap.update(pileupWeight2018)

pileupWeight = cms.EDProducer("CATPileupWeightProducer",
    #weightingMethod = cms.string("NVertex"), # Simple bin-by-bin correction of nVertex distribution. Non standard
    weightingMethod = cms.string("Standard"), # The Standard method in the CMSSW
    #weightingMethod = cms.string("RedoWeight"), # this is to be used re-reweight on CATTuple
    pileupMC = cms.vdouble(),
    pileupRD = cms.vdouble(),
    pileupUp = cms.vdouble(),
    pileupDn = cms.vdouble(),
    simpleWeights = cms.vdouble(),
    #pileupInfo = cms.InputTag("addPileupInfo"), # For the AOD and MiniAODv1
    pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    nTrueIntr = cms.InputTag("pileupWeight", "nTrueInteraction", "CAT"),
    #LuminositySectionsBlockRange = LumiList('LumiMask/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()
)


