import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patSequences_cff import * 
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff import * 




#change to 8E29 menu
#patTrigger.processName      = 'HLT8E29'
#patTriggerEvent.processName = 'HLT8E29'

#Keep the matchingInfo for two paths: one physics and one backup
#electronTriggerMatchHLTEle15SWL1R = electronTriggerMatchHLTEle15LWL1R.clone()
#electronTriggerMatchHLTEle15SWL1R.maxDeltaR = 0.3
#electronTriggerMatchHLTEle15SWL1R.src = 'selectedPatElectrons'
#electronTriggerMatchHLTEle15SWL1R.pathNames = 'HLT_Ele15_SW_L1R',#backup
#electronTriggerMatchHLTEle10SWL1R = electronTriggerMatchHLTEle15SWL1R.clone()
#electronTriggerMatchHLTEle10SWL1R.pathNames = 'HLT_Ele10_SW_L1R',#physics

# keep matching only for electron
#patTriggerMatcher = cms.Sequence(electronTriggerMatchHLTEle10SWL1R+electronTriggerMatchHLTEle15SWL1R+electronTriggerMatchHLTDoubleEle5SWL1R)

#change the matching sources in embedder
#cleanPatElectronsTriggerMatch.matches = 'electronTriggerMatchHLTEle10SWL1R','electronTriggerMatchHLTEle15SWL1R'
#cleanPatElectronsTriggerMatch.src = 'selectedPatElectrons'

# keep Embedding only for electron
#patTriggerMatchEmbedder = cms.Sequence(cleanPatElectronsTriggerMatch)

# change sources. Otherwis it looks for muonMatches, etc.
#patTriggerEvent.patTriggerMatches = 'electronTriggerMatchHLTEle10SWL1R','electronTriggerMatchHLTEle15SWL1R'

#patDefaultWithTrigger = cms.Sequence(patDefaultSequence* patTrigger * patTriggerMatcher * patTriggerMatchEmbedder * patTriggerEvent)
#patAddTrigger = cms.Sequence(patTrigger * patTriggerMatcher * patTriggerMatchEmbedder * patTriggerEvent)

patDefaultWithTrigger = cms.Sequence(patDefaultSequence)
#patAddTrigger = cms.Sequence(patTrigger * patTriggerEvent)


