import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import listModules, applyPostfix
from CommonTools.ParticleFlow.pfParticleSelection_cff import *
from CommonTools.ParticleFlow.Isolation.pfMuonIsolation_cff import *

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
