import FWCore.ParameterSet.Config as cms

def catTool(process, runOnMC=True, doSecVertex=True, useMiniAOD = True, bunchCrossing=25):
    catJetsSource = "slimmedJets"
    catGenJetsSource = "slimmedGenJets"
    catMETsSource = "slimmedMETs"
    catJetsPuppiSource = "slimmedJetsPuppi"
    catMETsPuppiSource = "slimmedMETsPuppi"
    catMuonsSource = "slimmedMuons"
    catElectronsSource = "slimmedElectrons"
    catPhotonsSource = "slimmedPhotons"
    catTausSource = "slimmedTaus"
    catVertexSource = "offlineSlimmedPrimaryVertices"
    catVertex = "catVertex"
    catMCsource = "prunedGenParticles"
    catBeamSpot = "offlineBeamSpot"
    catRho = "fixedGridRhoAll"
    ePidNames = cms.vstring()
    btagNames = cms.vstring("pfCombinedInclusiveSecondaryVertexV2BJetTags")

    process.nEventsTotal = cms.EDProducer("EventCountProducer")
    process.nEventsFiltered = cms.EDProducer("EventCountProducer")
    process.load("CATTools.CatProducer.catCandidates_cff")
#######################################################################
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    ## Hcal HBHE
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
    if bunchCrossing == 50:
        process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(True) 
        process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun1")

    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
        reverseDecision = cms.bool(False)
    )

    process.p = cms.Path(
        process.nEventsTotal*
        process.HBHENoiseFilterResultProducer* #produces HBHE bools
        process.ApplyBaselineHBHENoiseFilter*  #reject events based
        process.nEventsFiltered
    )
#######################################################################
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
# recompute the T1 PFMET
    jecUncertaintyFile = "CATTools/CatProducer/data/Summer15_25nsV5_DATA_UncertaintySources_AK4PFchs.txt"
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD( process, isData= not runOnMC, jecUncFile=jecUncertaintyFile)
# MET without HF
    process.noHFCands = cms.EDFilter("CandPtrSelector",
                                     src=cms.InputTag("packedPFCandidates"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )
    runMetCorAndUncFromMiniAOD( process, isData=not runOnMC, jecUncFile=jecUncertaintyFile, pfCandColl=cms.InputTag("noHFCands"), postfix="NoHF")
    process.catMETsNoHF = process.catMETs.clone()
    process.catMETsNoHF.src = cms.InputTag("slimmedMETsNoHF")
    ## no residuals currently available 
    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
        
    del process.slimmedMETs.t01Variation
    #del process.slimmedMETs.t1Uncertainties
    del process.slimmedMETs.tXYUncForRaw
    del process.slimmedMETs.tXYUncForT1
    del process.slimmedMETsNoHF.t01Variation
    #del process.slimmedMETsNoHF.t1Uncertainties
    del process.slimmedMETsNoHF.tXYUncForRaw
    del process.slimmedMETsNoHF.tXYUncForT1
    
#######################################################################    
# adding pfMVAMet https://twiki.cern.ch/twiki/bin/viewauth/CMS/MVAMet#Spring15_samples_with_25ns_50ns
# https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/RecoMET/METPUSubtraction/test/mvaMETOnMiniAOD_cfg.py
    process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
    process.ak4PFJets.doAreaFastjet = cms.bool(True)

    from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

    process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
    #process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons")
    process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
    process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")

    process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
    #process.puJetIdForPFMVAMEt.jets = cms.InputTag("ak4PFJets")
    process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")

    #process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    #process.ak4PFJetsForPFMVAMet = process.ak4PFJets.clone()
    #process.ak4PFJetsForPFMVAMet.src = cms.InputTag("packedPFCandidates")
    #This is temporary solution to avoid the circular dependenc error. Hope the recipe for miniAOD is available soon.(Tae Jeong) 
    #https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideUnscheduledExecution#Circular_Dependence_Errors
    #from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3
    #process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
    #process.pfMVAMEt.srcUncorrJets = cms.InputTag("ak4PFJetsForPFMVAMet")
    #process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
    #process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
    #process.pfMVAMEt.srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt")
    #process.calibratedAK4PFJetsForPFMVAMEt.src = cms.InputTag("ak4PFJetsForPFMVAMet")
    #process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
    #process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
    #process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")
    #from RecoJets.JetProducers.PileupJetIDParams_cfi import full_74x_chs
    #process.puJetIdForPFMVAMEt.algos = cms.VPSet(full_74x_chs)

    process.pfMVAMEt.inputFileNames = cms.PSet(
        U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_25NS_July2015.root'),
        DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_25NS_July2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_25NS_July2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_25NS_July2015.root')
    )
        
    process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
    process.patMETsPfMva = process.patMETs.clone()
    process.patMETsPfMva.addGenMET    = cms.bool(False)
    process.patMETsPfMva.metSource  = cms.InputTag("pfMVAMEt")
    process.patMETsPfMva.muonSource = cms.InputTag(catMuonsSource)
    process.catMETsPfMva = process.catMETs.clone()
    process.catMETsPfMva.src = cms.InputTag("patMETsPfMva")
#######################################################################
# redoing puppi from miniAOD as recommended
# https://twiki.cern.ch/twiki/bin/view/CMS/PUPPI
    process.load('CommonTools/PileupAlgos/Puppi_cff')
    ## e.g. to run on miniAOD
    process.puppi.candName = cms.InputTag('packedPFCandidates')
    process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

    process.load('CommonTools/PileupAlgos/Puppi_cff')
    process.puppi.candName = cms.InputTag('packedPFCandidates')
    process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
    # remaking puppi jets
    from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
    jetToolbox( process, 'ak4', 'ak4JetSubs', 'out', PUMethod='Puppi', miniAOD = useMiniAOD, runOnMC = True,#due to bug in jetToolbox
                JETCorrPayload='AK4PFPuppi', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] )#bug-JETCorrLevels overwritten in jetToolbox

    process.selectedPatJetsAK4PFPuppi.cut = cms.string("pt > 20")
    process.patJetGenJetMatchAK4PFPuppi.matched = cms.InputTag("slimmedGenJets")
    process.patJetsAK4PFPuppi.embedGenPartonMatch = cms.bool(False)
    catJetsPuppiSource = "selectedPatJetsAK4PFPuppi"
    # remaking puppi met
    from RecoMET.METProducers.PFMET_cfi import pfMet
    process.pfMetPuppi = pfMet.clone();
    process.pfMetPuppi.src = cms.InputTag('puppi')
    process.patPfMetPuppi = process.patMETs.clone()
    process.patPfMetPuppi.addGenMET    = cms.bool(False)
    process.patPfMetPuppi.metSource  = cms.InputTag("pfMetPuppi")
    process.patPfMetPuppi.muonSource = cms.InputTag(catMuonsSource)
    catMETsPuppiSource = "patPfMetPuppi"

    # for puppi isolation
    ## process.packedPFCandidatesWoMuon  = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV>=2 && abs(pdgId)!=13 " ) )
    ## process.particleFlowNoMuonPUPPI.candName         = 'packedPFCandidatesWoMuon'
    ## process.particleFlowNoMuonPUPPI.vertexName       = 'offlineSlimmedPrimaryVertices'
    
#######################################################################    
# getting jec from file for jec on the fly from db file
# currently only for mc
    useJECfile = True
    era = "Summer15_25nsV5"
    if runOnMC:
        era = era+"_MC"
    else:
        era = era+"_DATA"
        
    if useJECfile:
        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
            connect = cms.string('sqlite_fip:CATTools/CatProducer/data/'+era+'.db'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                    label= cms.untracked.string("AK4PF")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                    label= cms.untracked.string("AK4PFchs")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFPuppi"),
                    label= cms.untracked.string("AK4PFPuppi")),
            )
        )
        process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
        print "JEC based on", process.jec.connect
    
## applying new jec on the fly
    if useMiniAOD:
        process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
        process.patJetCorrFactors.primaryVertices = cms.InputTag(catVertexSource)
        catJetsSource = "patJetsUpdated"
        ### updating puppi jet jec
        process.patJetPuppiCorrFactorsUpdated = process.patJetCorrFactorsUpdated.clone(
        payload = cms.string('AK4PFPuppi'),
            src = cms.InputTag(catJetsPuppiSource),
        )
        process.patJetsPuppiUpdated = process.patJetsUpdated.clone(
            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetPuppiCorrFactorsUpdated")),
            jetSource = cms.InputTag(catJetsPuppiSource),
        )
        catJetsPuppiSource = "patJetsPuppiUpdated"

#######################################################################
## for egamma pid temp
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_74X
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection
    my_id_modules = [
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                     #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                     ]
    if useMiniAOD:
        dataFormat = DataFormat.MiniAOD
    else :
        dataFormat = DataFormat.AOD

    switchOnVIDElectronIdProducer(process, dataFormat)

    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
    process.catElectrons.electronIDSources = cms.PSet(
        eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
        eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
        eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
        mvaNonTrigMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
        mvaNonTrigTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
        #mvaTrigMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
        #mvaTrigTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
    )

#######################################################################    
    if runOnMC:## Load MC dependent producers
        ## FIX ME - pile up and pdf weight
        process.load("CATTools.CatProducer.genWeight_cff")
        process.load("CATTools.CatProducer.pileupWeight_cff")
        process.pileupWeight.vertex = cms.InputTag(catVertexSource)

        process.catMuons.shiftedEnDownSrc = cms.InputTag("shiftedPatMuonEnDown")
        process.catMuons.shiftedEnUpSrc = cms.InputTag("shiftedPatMuonEnUp")
        process.catElectrons.shiftedEnDownSrc = cms.InputTag("shiftedPatElectronEnDown")
        process.catElectrons.shiftedEnUpSrc = cms.InputTag("shiftedPatElectronEnUp")

    if doSecVertex:
        from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
        setattr(process, "TransientTrackBuilderESProducer", TransientTrackBuilderESProducer)
        #process.makeCatCandidates += process.catSecVertexs
                
    process.catJets.src = cms.InputTag(catJetsSource)
    process.catJets.btagNames = btagNames
    process.catJets.payloadName = cms.string("AK4PFchs")
    process.catTaus.src = cms.InputTag(catTausSource)
    process.catMuons.src = cms.InputTag(catMuonsSource)
    process.catMuons.mcLabel = cms.InputTag(catMCsource)
    process.catMuons.vertexLabel = cms.InputTag(catVertex)
    process.catMuons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.src = cms.InputTag(catElectronsSource)
    process.catElectrons.ePidNames = ePidNames
    process.catElectrons.vertexLabel = cms.InputTag(catVertex)
    process.catElectrons.mcLabel = cms.InputTag(catMCsource)
    process.catElectrons.beamLineSrc = cms.InputTag(catBeamSpot)
    process.catElectrons.rhoLabel = cms.InputTag(catRho)
    process.catPhotons.src = cms.InputTag(catPhotonsSource)
    process.catMETs.src = cms.InputTag(catMETsSource)
    process.catSecVertexs.muonSrc = cms.InputTag(catMuonsSource)
    process.catSecVertexs.elecSrc = cms.InputTag(catElectronsSource)
    process.catSecVertexs.vertexLabel = cms.InputTag(catVertexSource)

    process.catJetsPuppi.src = cms.InputTag(catJetsPuppiSource)
    process.catJetsPuppi.btagNames = btagNames
    process.catJetsPuppi.payloadName = cms.string("AK4PFchs") #temp for now
    process.catMETsPuppi.src = cms.InputTag(catMETsPuppiSource)
    process.catVertex.vertexLabel = cms.InputTag(catVertexSource)
    
