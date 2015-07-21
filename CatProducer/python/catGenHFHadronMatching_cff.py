import FWCore.ParameterSet.Config as cms

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
        from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJets
        process.genParticlesForJets = genParticlesForJets.clone()
   
    # Supplies PDG ID to real name resolution of MC particles
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    
    # Producing own jets for testing purposes
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    process.ak4GenJetsCustom = ak4GenJets.clone(
        src = genJetInputParticleCollection,
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
    )
   
    # Ghost particle collection used for Hadron-Jet association 
    # MUST use proper input particle collection
    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
    process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
        particles = genParticleCollection
    )
   
    # Input particle collection for matching to gen jets (partons + leptons) 
    # MUST use use proper input jet collection: the jets to which hadrons should be associated
    # rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
    # More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition
    from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import genJetFlavourPlusLeptonInfos
    process.genJetFlavourPlusLeptonInfos = genJetFlavourPlusLeptonInfos.clone(
        jets = genJetCollection,
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
    )
   
   
    # Plugin for analysing B hadrons
    # MUST use the same particle collection as in selectedHadronsAndPartons
    from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import matchGenBHadron
    process.matchGenBHadron = matchGenBHadron.clone(
        genParticles = genParticleCollection
    )
   
    # Plugin for analysing C hadrons
    # MUST use the same particle collection as in selectedHadronsAndPartons
    from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import matchGenCHadron
    process.matchGenCHadron = matchGenCHadron.clone(
        genParticles = genParticleCollection
    )
    
