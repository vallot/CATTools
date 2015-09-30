import FWCore.ParameterSet.Config as cms

#ref. https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/PhysicsTools/JetMCAlgos/test/matchGenHFHadrons.py#L60
def genHFTool(process, useMiniAOD = True):
    # Setting input particle collections to be used by the tools
    genParticleCollection = ''
    genJetInputParticleCollection = ''
    genJetCollection = 'ak4GenJetsCustom'
    if useMiniAOD:
        genParticleCollection = 'prunedGenParticles'
        genJetInputParticleCollection = 'packedGenParticles'
    else:
        genParticleCollection = 'genParticles'
        genJetInputParticleCollection = 'genParticlesForJets'
        ## producing a subset of genParticles to be used for jet clustering in AOD
        process.load("RecoJets.Configuration.GenJetParticles_cff")
   
    # Supplies PDG ID to real name resolution of MC particles
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    
    # Producing own jets for testing purposes
    process.load("RecoJets.JetProducers.ak4GenJets_cfi")
    process.ak4GenJetsCustom = process.ak4GenJets.clone(
        src = cms.InputTag(genJetInputParticleCollection),
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
    )
   
    # Ghost particle collection used for Hadron-Jet association 
    # MUST use proper input particle collection
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.selectedHadronsAndPartons.particles = genParticleCollection
   
    # Input particle collection for matching to gen jets (partons + leptons) 
    # MUST use use proper input jet collection: the jets to which hadrons should be associated
    # rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
    # More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition
    process.load("PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff")
    process.genJetFlavourPlusLeptonInfos.jets = genJetCollection
    process.genJetFlavourPlusLeptonInfos.rParam = cms.double(0.4)
    process.genJetFlavourPlusLeptonInfos.jetAlgorithm = cms.string("AntiKt")
   
    # Plugin for analysing B hadrons
    # MUST use the same particle collection as in selectedHadronsAndPartons
    process.load("PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff")
    process.matchGenBHadron.genParticles = genParticleCollection
   
    # Plugin for analysing C hadrons
    # MUST use the same particle collection as in selectedHadronsAndPartons
    process.load("PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff")
    process.matchGenCHadron.genParticles = genParticleCollection

    process.load("CATTools.CatProducer.mcTruthTop.GenTtbarCategorizer_cfi")
    process.GenTtbarCategories = process.categorizeGenTtbar.clone(
       genJets = cms.InputTag(genJetCollection)
    )
