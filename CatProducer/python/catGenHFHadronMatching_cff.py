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
       genJets = cms.InputTag(genJetCollection)
    )
   
