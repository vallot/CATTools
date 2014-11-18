import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import listModules, applyPostfix
from CommonTools.ParticleFlow.pfParticleSelection_cff import *
from CommonTools.ParticleFlow.Isolation.pfMuonIsolation_cff import *

def miniAOD(process):
  process.catMuons.src = "slimmedMuons" 
  process.catMuons.mcLabel = "prunedGenParticles"
  process.catMuons.vertexLabel = "offlineSlimmedPrimaryVertices" 
  process.catElectrons.src = "slimmedElectrons"
  process.catElectrons.mcLabel = "prunedGenParticles"
  process.catElectrons.vertexLabel = "offlineSlimmedPrimaryVertices"
  process.catPhotons.src = "slimmedPhotons"
  process.catMETs.src = "slimmedMETs"
  process.catJets.src = "slimmedJets"
  process.catGenJets.src = "slimmedGenJets" 
  process.catMCParticles.src = "prunedGenParticles"
  process.catGenTops.genJetLabel = "slimmedGenJets"
  process.catGenTops.mcParticleLabel = "prunedGenParticles" 

#muons with weighted isolation method
def addMuonWeighted(process):
  ## for muon isolation study with weighted method
  print "4508 PR needs to be merged"
  print "git cms-merge-topic 4508"
  print "or use the release above CMSSW_7_0_7"

  ### add different cone size
  #let's use userIsolation function to use different cone size for the time being
  #userIsolation("pat::User1Iso") for chargedHadronIso()
  #userIsolation("pat::User2Iso") for neutralHadronIso()
  #userIsolation("pat::User3Iso") for photonIso()
  #userIsolation("pat::User4Iso") for puChargedHadronIso()
  #userIsolation("pat::User5Iso") for particleIso()
  #for muon
  process.patMuons.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03","muPFIsoValueNeutral03","muPFIsoValueGamma03","muPFIsoValuePU03","muPFIsoValueChargedAll03")
  #for electron
  process.patElectrons.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFId","elPFIsoValueNeutral03PFId","elPFIsoValueGamma03PFId","elPFIsoValuePU03PFId","elPFIsoValueChargedAll03PFId")
  process.patElectrons.isolationValuesNoPFId.user = cms.VInputTag("elPFIsoValueCharged03NoPFId","elPFIsoValueNeutral03NoPFId","elPFIsoValueGamma03NoPFId","elPFIsoValuePU03NoPFId","elPFIsoValueChargedAll03NoPFId")
  ###

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
      pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralWeighted" ),
      pfPhotons = cms.InputTag("muPFIsoDepositGammaWeighted" ),
      pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPUWeighted" ),
      pfChargedAll = cms.InputTag("muPFIsoDepositChargedAllWeighted" ),
      )

  process.patMuonsWeighted.isolationValues = cms.PSet(
      pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04Weighted"),
      pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04Weighted" ),
      pfPhotons = cms.InputTag("muPFIsoValueGamma04Weighted" ),
      pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04Weighted" ),
      pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04Weighted"),
      user = cms.VInputTag("muPFIsoValueCharged03Weighted","muPFIsoValueNeutral03Weighted","muPFIsoValueGamma03Weighted","muPFIsoValuePU03Weighted","muPFIsoValueChargedAll03Weighted"),
      )


def useRecoMuon(process, postfix, dR = "03"):
 
  sourceMuons = 'muons' #take it from reco collection

  module = applyPostfix(process,"patMuons",postfix)

  module.pfMuonSource = sourceMuons
  module.useParticleFlow = False

  muPFIsoDepositCharged.src = sourceMuons
  muPFIsoDepositChargedAll.src = sourceMuons
  muPFIsoDepositNeutral.src = sourceMuons
  muPFIsoDepositGamma.src = sourceMuons
  muPFIsoDepositPU.src = sourceMuons

  module.isoDeposits = cms.PSet(
      pfChargedHadrons = cms.InputTag("muPFIsoDepositCharged" ),
      pfChargedAll = cms.InputTag("muPFIsoDepositChargedAll" ),
      pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPU" ),
      pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutral" ),
      pfPhotons = cms.InputTag("muPFIsoDepositGamma" ),
      )

  module.isolationValues = cms.PSet(
      pfChargedHadrons = cms.InputTag("muPFIsoValueCharged"+dR),
      pfChargedAll = cms.InputTag("muPFIsoValueChargedAll"+dR),
      pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU"+dR ),
      pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral"+dR ),
      pfPhotons = cms.InputTag("muPFIsoValueGamma"+dR ),
      )

  applyPostfix(process,"muonMatch",postfix).src = module.pfMuonSource
  getattr(process,'patDefaultSequence'+postfix).replace( getattr(process,"patMuons"+postfix),
                                                   pfParticleSelectionSequence+
                                                   pfMuonIsolationSequence+
                                                   getattr(process,"patMuons"+postfix) )
