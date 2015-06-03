import FWCore.ParameterSet.Config as cms

def catTool(process, runOnMC=True, doSecVertex=True, useMiniAOD = True):
    catJetsSource = "slimmedJets"
    catGenJetsSource = "slimmedGenJets"
    catMuonsSource = "slimmedMuons"
    catElectronsSource = "slimmedElectrons"
    catPhotonsSource = "slimmedPhotons"
    catTausSource = "slimmedTaus"
    catMETsSource = "slimmedMETs"
    catVertexSource = "offlineSlimmedPrimaryVertices"
    catMCsource = "prunedGenParticles"
    catBeamSpot = "offlineBeamSpot"
    catRho = "fixedGridRhoAll"

    process.load("CATTools.CatProducer.catCandidates_cff")        
    process.load("CATTools.CatProducer.recoEventInfo_cfi")

    if runOnMC:## Load MC dependent producers
        process.load("CATTools.CatProducer.pdfWeight_cff")
        process.load("CATTools.CatProducer.pileupWeight_cff")
        process.load("CATTools.CatProducer.pseudoTop_cfi")
        if not useMiniAOD:
            process.load("CATTools.CatProducer.genTopProducer_cfi")

    process.catJets.src = cms.InputTag(catJetsSource)
    process.catJets.genJetMatch = cms.InputTag("patJetGenJetMatch")
    process.catTaus.src = cms.InputTag(catTausSource)
    process.catTaus.genJetMatch = cms.InputTag("tauGenJetMatch")
    process.catMuons.src = cms.InputTag(catMuonsSource)
    process.catMuons.mcLabel = cms.InputTag(catMCsource)
    process.catMuons.vertexLabel = cms.InputTag(catVertexSource)
    process.catMuons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.src = cms.InputTag(catElectronsSource)
    process.catElectrons.vertexLabel = cms.InputTag(catVertexSource)
    process.catElectrons.mcLabel = cms.InputTag(catMCsource)
    process.catElectrons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.rhoLabel = cms.InputTag(catRho)
    process.catPhotons.src = cms.InputTag(catPhotonsSource)
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catSecVertexs.muonSrc = cms.InputTag(catMuonsSource)
    process.catSecVertexs.elecSrc = cms.InputTag(catElectronsSource)
    process.catSecVertexs.vertexLabel = cms.InputTag(catVertexSource)

    if runOnMC:
        ## using MEtUncertainties to get lepton shifts
        ## Need to update - Jet calculation was old, would most likely be the same for leptons
        from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
        runType1PFMEtUncertainties(process,
                                    addToPatDefaultSequence=False,
                                    jetCollection=catJetsSource,
                                    electronCollection=catElectronsSource,
                                    photonCollection=catPhotonsSource,
                                    muonCollection=catMuonsSource,
                                    tauCollection=catTausSource,
                                    makeType1p2corrPFMEt=True,
                                    outputModule=None)
        
        process.catMuons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedMuonsEnDown")
        process.catMuons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedMuonsEnUp")
        process.catElectrons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedElectronsEnDown")
        process.catElectrons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedElectronsEnUp")

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        #process.makeCatCandidates += process.catSecVertexs

    ## for egamma pid temp 
    ## https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_74X
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection
    if not useMiniAOD :
        dataFormat = DataFormat.AOD
    else :
        dataFormat = DataFormat.MiniAOD
            
    switchOnVIDElectronIdProducer(process, dataFormat)
    
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']
    
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
    if useMiniAOD:
        ## applying new jec on the fly
        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
            connect = cms.string('sqlite_file:../data/PHYS14_V4_MC.db'),
            toGet = cms.VPSet(
                cms.PSet(record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_PHYS14_V4_MC_AK4PF"),
                label= cms.untracked.string("AK4PF"))
            ))
        process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

        #process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
        #process.catJets.src = cms.InputTag("patJetsUpdated")

        from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

        process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
    
        process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)
        process.ak4PFJetsCHS = ak4PFJets.clone(src = 'chs', doAreaFastjet = True)
        process.ak4GenJets = ak4GenJets.clone(src = 'packedGenParticles')

        process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
        from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
        addJetCollection(
            process,
            postfix   = "",
            labelName = 'AK4PFCHS',
            jetSource = cms.InputTag('ak4PFJetsCHS'),
            pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
            pfCandidates = cms.InputTag('packedPFCandidates'),
            svSource = cms.InputTag('slimmedSecondaryVertices'),
            jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
            btagDiscriminators = [ 'pfCombinedSecondaryVertexBJetTags', 'pfCombinedInclusiveSecondaryVertexV2BJetTags' ],
            genJetCollection=cms.InputTag('ak4GenJets'),
            genParticles=cms.InputTag('prunedGenParticles')
        )


        process.catJets.btagNames = cms.vstring()
        process.catJets.src = cms.InputTag("patJetsAK4PFCHS")

        process.catElectrons.electronIDSources = cms.PSet(
            eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
            eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
            eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
            eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
            eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
        )
