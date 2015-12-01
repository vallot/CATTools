import FWCore.ParameterSet.Config as cms

def catTool(process, runOnMC=True, useMiniAOD=True):
    bunchCrossing=25
    globaltag_run2_50ns = ["MCRUN2_74_V9A", "74X_mcRun2_startup_v2", "74X_mcRun2_asymptotic50ns_v0", "74X_dataRun2_v2"]
    for i in globaltag_run2_50ns:
        if i == process.GlobalTag.globaltag:
            bunchCrossing=50

    if runOnMC and bunchCrossing == 50:
        from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
        process.pileupWeight.pileupMC = pileupWeightMap["Startup2015_50ns"]
        process.pileupWeight.pileupRD = pileupWeightMap["Run2015_50ns"]
        process.pileupWeight.pileupUp = pileupWeightMap["Run2015_50nsUp"]
        process.pileupWeight.pileupDn = pileupWeightMap["Run2015_50nsDn"]

    lumiJSON = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON'
    if runOnMC:
        from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
        process.pileupWeight.pileupMC = pileupWeightMap["Startup2015_25ns"]
        process.pileupWeight.pileupRD = pileupWeightMap["%s"%lumiJSON]
        process.pileupWeight.pileupUp = pileupWeightMap["%s_Up"%lumiJSON]
        process.pileupWeight.pileupDn = pileupWeightMap["%s_Dn"%lumiJSON]
        process.pileupWeightSilver = process.pileupWeight.clone()
        process.pileupWeightSilver.pileupRD = pileupWeightMap["%s_Silver"%lumiJSON]
        process.pileupWeightSilver.pileupUp = pileupWeightMap["%s_Silver_Up"%lumiJSON]
        process.pileupWeightSilver.pileupDn = pileupWeightMap["%s_Silver_Dn"%lumiJSON]
    else:
        from FWCore.PythonUtilities.LumiList import LumiList
        process.lumiMask = cms.EDProducer("LumiMaskProducer",
            LumiSections = LumiList('../data/LumiMask/%s.txt'%lumiJSON).getVLuminosityBlockRange())
        process.lumiMaskSilver = cms.EDProducer("LumiMaskProducer",
            LumiSections = LumiList('../data/LumiMask/%s_Silver.txt'%lumiJSON).getVLuminosityBlockRange())
    
    useJECfile = True
    era = "Summer15_{}nsV6".format(bunchCrossing)
    
    jecUncertaintyFile = "CATTools/CatProducer/data/JEC/%s_DATA_UncertaintySources_AK4PFchs.txt"%era
    if runOnMC:
        era = era+"_MC"
    else:
        era = era+"_DATA"
    
    if useJECfile:
        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
            connect = cms.string('sqlite_fip:CATTools/CatProducer/data/JEC/%s.db'%era),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_%s_AK4PF"%era),
                    label= cms.untracked.string("AK4PF")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFchs"%era),
                    label= cms.untracked.string("AK4PFchs")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFPuppi"%era),
                    label= cms.untracked.string("AK4PFPuppi")),
            )
        )
        process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
        print "JEC based on", process.jec.connect
    
    if useMiniAOD: ## corrections when using miniAOD
        #######################################################################
        ## Hcal HBHE https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
        process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
        process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
        process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
        process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

        if bunchCrossing == 50:
            process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(True) 
            process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun1")

        process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
            inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
            reverseDecision = cms.bool(False))
        
        process.nEventsFiltered = cms.EDProducer("EventCountProducer")
    
        process.p += (process.HBHENoiseFilterResultProducer* #produces HBHE bools
                      process.ApplyBaselineHBHENoiseFilter*  #reject events based
                      process.nEventsFiltered)
        #######################################################################
        # adding puppi https://twiki.cern.ch/twiki/bin/view/CMS/PUPPI        
        #process.catJetsPuppi.src = cms.InputTag("slimmedJetsPuppi")
        #process.catMETsPuppi.src = cms.InputTag("slimmedMETsPuppi")
        # for puppi isolation
        ## process.packedPFCandidatesWoMuon  = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV>=2 && abs(pdgId)!=13 " ) )
        ## process.particleFlowNoMuonPUPPI.candName         = 'packedPFCandidatesWoMuon'
        ## process.particleFlowNoMuonPUPPI.vertexName       = 'offlineSlimmedPrimaryVertices'
        
        #######################################################################
        ## applying new jec on the fly
        process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
        process.patJetCorrFactors.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
        process.catJets.src = cms.InputTag("patJetsUpdated")

        ### updating puppi jet jec
        process.patJetPuppiCorrFactorsUpdated = process.patJetCorrFactorsUpdated.clone(
            src = process.catJetsPuppi.src,
            payload = cms.string('AK4PFPuppi'),
            levels = cms.vstring('L2Relative','L3Absolute'),
            useRho = cms.bool(False))
        
        process.patJetsPuppiUpdated = process.patJetsUpdated.clone(
            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetPuppiCorrFactorsUpdated")),
            jetSource = process.catJetsPuppi.src )
        
        process.catJetsPuppi.src = cms.InputTag("patJetsPuppiUpdated")
        process.catJetsPuppi.setGenParticle = cms.bool(False)
        #######################################################################
        # MET corrections from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
        runMetCorAndUncFromMiniAOD( process, isData= not runOnMC, jecUncFile=jecUncertaintyFile, jetColl= process.catJets.src)
        process.catMETs.src = cms.InputTag("slimmedMETs","","CAT")
        del process.slimmedMETs.caloMET

        ## redoing noHF met due to new correction
        process.noHFCands = cms.EDFilter("CandPtrSelector",src=cms.InputTag("packedPFCandidates"),
                                         cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0"))
        runMetCorAndUncFromMiniAOD(process,isData=not runOnMC,pfCandColl=cms.InputTag("noHFCands"),postfix="NoHF",
                                   jecUncFile=jecUncertaintyFile, jetColl= process.catJets.src)

        process.catMETsNoHF = process.catMETs.clone(src = cms.InputTag("slimmedMETsNoHF","","CAT"))
        del process.slimmedMETsNoHF.caloMET
        #######################################################################
        ## for egamma pid https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_74X
        from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection
        electron_ids = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                        'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']
        switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
        for idmod in electron_ids:
            setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

        process.catElectrons.electronIDSources = cms.PSet(
            cutBasedElectronID_Spring15_25ns_V1_standalone_loose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
            cutBasedElectronID_Spring15_25ns_V1_standalone_medium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
            cutBasedElectronID_Spring15_25ns_V1_standalone_tight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
            cutBasedElectronID_Spring15_25ns_V1_standalone_veto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
            heepElectronID_HEEPV60 = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
            mvaEleID_Spring15_25ns_nonTrig_V1_wp80 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
            mvaEleID_Spring15_25ns_nonTrig_V1_wp90 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
            mvaEleID_Spring15_25ns_Trig_V1_wp80 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
            mvaEleID_Spring15_25ns_Trig_V1_wp90 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
        )
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

    process.pfMVAMEt.inputFileNames = cms.PSet(
        U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_{}NS_July2015.root'.format(bunchCrossing)),
        DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_{}NS_July2015.root'.format(bunchCrossing)),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_{}NS_July2015.root'.format(bunchCrossing)),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_{}NS_July2015.root'.format(bunchCrossing))
    )
        
    process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
    process.patMETsPfMva = process.patMETs.clone(addGenMET = cms.bool(False), metSource  = cms.InputTag("pfMVAMEt"))
    process.catMETsPfMva = process.catMETs.clone(src = cms.InputTag("patMETsPfMva"))
    process.catMETsPfMva.setUnclusteredEn = cms.bool(False)
    process.catMETsPfMva.setJetMETSyst = cms.bool(False)
    
