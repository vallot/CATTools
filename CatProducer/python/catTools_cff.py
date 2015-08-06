import FWCore.ParameterSet.Config as cms

def catTool(process, runOnMC=True, doSecVertex=True, useMiniAOD = True):
    catJetsSource = "slimmedJets"
    catJetsPuppiSource = "slimmedJetsPuppi"
    catGenJetsSource = "slimmedGenJets"
    catMETsSource = "slimmedMETs"
    catMETsPuppiSource = "slimmedMETsPuppi"
    catMuonsSource = "slimmedMuons"
    catElectronsSource = "slimmedElectrons"
    catPhotonsSource = "slimmedPhotons"
    catTausSource = "slimmedTaus"
    catVertexSource = "offlineSlimmedPrimaryVertices"
    catMCsource = "prunedGenParticles"
    catBeamSpot = "offlineBeamSpot"
    catRho = "fixedGridRhoAll"
    btagNames = cms.vstring("pfCombinedInclusiveSecondaryVertexV2BJetTags")
    ePidNames = cms.vstring()

    process.nEventsTotal = cms.EDProducer("EventCountProducer")
    #process.p = cms.Path(process.nEventsTotal)
    process.load("CATTools.CatProducer.catCandidates_cff")        
#######################################################################    
# adding pfMVAMet
    process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    process.ak4PFJetsForPFMVAMet = process.ak4PFJets.clone()
    process.ak4PFJetsForPFMVAMet.src = cms.InputTag("packedPFCandidates")
    #This is temporary solution to avoid the circular dependenc error. Hope the recipe for miniAOD is available soon.(Tae Jeong) 
    #https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideUnscheduledExecution#Circular_Dependence_Errors
    from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3
    process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
    process.pfMVAMEt.srcUncorrJets = cms.InputTag("ak4PFJetsForPFMVAMet")
    process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
    process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
    process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")
    process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
    process.patMETsPfMva = process.patMETs.clone()
    process.patMETsPfMva.addGenMET    = cms.bool(False)
    process.patMETsPfMva.metSource  = cms.InputTag("pfMVAMEt")
    process.patMETsPfMva.muonSource = cms.InputTag(catMuonsSource)
    process.catMETsPfMva = process.catMETs.clone()
    process.catMETsPfMva.src = cms.InputTag("patMETsPfMva")
#######################################################################
#######################################################################    
# getting jec from file for jec on the fly from db file
# currently only for mc
    if runOnMC:
        #era = "PHYS14_V4_MC"
        era = "Summer15_50nsV2_MC"
        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
            connect = cms.string('sqlite_fip:CATTools/CatProducer/data/'+era+'.db'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                    label= cms.untracked.string("AK4PFchs")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PUPPI"),
                    label= cms.untracked.string("AK4PUPPI")),
            ))
        process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
## applying new jec on the fly
        if useMiniAOD:
            process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
            catJetsSource = "patJetsUpdated"
            ### updating puppi jet jec
            process.patJetPuppiCorrFactorsUpdated = process.patJetCorrFactorsUpdated.clone(
                payload = cms.string('AK4PUPPI'),
                src = cms.InputTag(catJetsPuppiSource),
            )
            process.patJetsPuppiUpdated = process.patJetsUpdated.clone(
                jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetPuppiCorrFactorsUpdated")),
                jetSource = cms.InputTag(catJetsPuppiSource),
            )
            catJetsPuppiSource = "patJetsPuppiUpdated"

#######################################################################
#######################################################################    
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
        process.catElectrons.electronIDSources = cms.PSet(
            eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
            eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
            eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
            eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
            eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
        )
#######################################################################    
    if runOnMC:## Load MC dependent producers
        ## FIX ME - pile up and pdf weight
        process.load("CATTools.CatProducer.pdfWeight_cff")
        process.load("CATTools.CatProducer.pileupWeight_cff")
        process.pileupWeight.vertex = cms.InputTag(catVertexSource)

        if not useMiniAOD:
            process.load("CATTools.CatProducer.genTopProducer_cfi")
            
        ## FIX ME very out of date!
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
                                    outputModule=None,
                                    jecUncertaintyFile=None,
                                    )
        
        process.catMuons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedMuonsEnDown")
        process.catMuons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedMuonsEnUp")
        process.catElectrons.shiftedEnDownSrc = cms.InputTag("shiftedSlimmedElectronsEnDown")
        process.catElectrons.shiftedEnUpSrc = cms.InputTag("shiftedSlimmedElectronsEnUp")

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        #process.makeCatCandidates += process.catSecVertexs
                
    process.catJets.src = cms.InputTag(catJetsSource)
    process.catJets.genJetMatch = cms.InputTag("patJetGenJetMatch")
    process.catJets.btagNames = btagNames
    process.catTaus.src = cms.InputTag(catTausSource)
    process.catTaus.genJetMatch = cms.InputTag("tauGenJetMatch")
    process.catMuons.src = cms.InputTag(catMuonsSource)
    process.catMuons.mcLabel = cms.InputTag(catMCsource)
    process.catMuons.vertexLabel = cms.InputTag(catVertexSource)
    process.catMuons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.src = cms.InputTag(catElectronsSource)
    process.catElectrons.ePidNames = ePidNames
    process.catElectrons.vertexLabel = cms.InputTag(catVertexSource)
    process.catElectrons.mcLabel = cms.InputTag(catMCsource)
    process.catElectrons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.rhoLabel = cms.InputTag(catRho)
    process.catPhotons.src = cms.InputTag(catPhotonsSource)
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catSecVertexs.muonSrc = cms.InputTag(catMuonsSource)
    process.catSecVertexs.elecSrc = cms.InputTag(catElectronsSource)
    process.catSecVertexs.vertexLabel = cms.InputTag(catVertexSource)

    process.catJetsPuppi.src = cms.InputTag(catJetsPuppiSource)
    process.catJetsPuppi.genJetMatch = cms.InputTag("patJetGenJetMatch")
    process.catJetsPuppi.btagNames = btagNames
    process.catMETsPuppi.src = cms.InputTag(catMETsPuppiSource)
    process.catVertex.vertexLabel = cms.InputTag(catVertexSource)
