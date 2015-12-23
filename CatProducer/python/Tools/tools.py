import FWCore.ParameterSet.Config as cms

#ref. https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/PhysicsTools/JetMCAlgos/test/matchGenHFHadrons.py#L60
def genHFTool(process, useMiniAOD = True):
    # Setting input particle collections to be used by the tools
    genParticleCollection = ''
    genJetCollection = 'ak4GenJetsCustom'
    if useMiniAOD:
        genParticleCollection = 'prunedGenParticles'
        genJetCollection = 'slimmedGenJets'
    else:
        genParticleCollection = 'genParticles'
        ## producing a subset of genParticles to be used for jet reclustering
        from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
        process.genParticlesForJetsCustom = genParticlesForJetsNoNu.clone(
            src = genParticleCollection
        )
        # Producing own jets for testing purposes
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        process.ak4GenJetsCustom = ak4GenJets.clone(
            src = 'genParticlesForJetsCustom',
            rParam = cms.double(0.4),
            jetAlgorithm = cms.string("AntiKt")
        )
    
    # Supplies PDG ID to real name resolution of MC particles
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

    # Ghost particle collection used for Hadron-Jet association 
    # MUST use proper input particle collection
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.selectedHadronsAndPartons.particles = genParticleCollection
    
    ## Input particle collection for matching to gen jets (partons + leptons) 
    # MUST use use proper input jet collection: the jets to which hadrons should be associated
    # rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
    # More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition
    from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
    process.genJetFlavourInfos = ak4JetFlavourInfos.clone( jets = genJetCollection )

    #for cmssw_7_6_X, not ready for 74x
    #process.load("PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff")
    #added the 3-lines instead of GenHFHadronMatcher_cff
    from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cfi import matchGenHFHadron
    process.matchGenBHadron=matchGenHFHadron.clone( flavour = 5)
    process.matchGenCHadron=matchGenHFHadron.clone( flavour = 4)   
    ## Plugin for analysing B hadrons
    # MUST use the same particle collection as in selectedHadronsAndPartons
    process.matchGenBHadron.genParticles    = genParticleCollection
    process.matchGenBHadron.jetFlavourInfos = "genJetFlavourInfos"
    ## Plugin for analysing C hadrons
    # MUST use the same particle collection as in selectedHadronsAndPartons
    process.matchGenCHadron.genParticles = genParticleCollection
    process.matchGenCHadron.jetFlavourInfos = "genJetFlavourInfos"
     
    process.load("TopQuarkAnalysis.TopTools.GenTtbarCategorizer_cfi")
    process.GenTtbarCategories = process.categorizeGenTtbar.clone(
       genJets = cms.InputTag(genJetCollection),
       genJetPtMin = cms.double(20.)
    )
    process.GenTtbarCategories30 = process.categorizeGenTtbar.clone(
       genJets = cms.InputTag(genJetCollection),
       genJetPtMin = cms.double(30.)
    )
    process.GenTtbarCategories40 = process.categorizeGenTtbar.clone(
       genJets = cms.InputTag(genJetCollection),
       genJetPtMin = cms.double(40.)
    )

    process.catGenTops.genJetLabel = genJetCollection
    process.catGenTops.mcParticleLabel = genParticleCollection


from PhysicsTools.PatAlgos.tools.helpers import listModules, applyPostfix
from CommonTools.ParticleFlow.pfParticleSelection_cff import *
from CommonTools.ParticleFlow.Isolation.pfMuonIsolation_cff import *

def miniAODWithoutGen(process) :
  process.catMuons.src = "slimmedMuons" 
  process.catMuons.mcLabel = "prunedGenParticles"
  process.catMuons.vertexLabel = "offlineSlimmedPrimaryVertices" 
  process.catElectrons.src = "slimmedElectrons"
  process.catElectrons.mcLabel = "prunedGenParticles"
  process.catElectrons.vertexLabel = "offlineSlimmedPrimaryVertices"
  process.catPhotons.src = "slimmedPhotons"
  process.catMETs.src = "slimmedMETs"
  process.catJets.src = "slimmedJets"


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
