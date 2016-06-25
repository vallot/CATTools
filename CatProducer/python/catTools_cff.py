import FWCore.ParameterSet.Config as cms
import catDefinitions_cfi as cat
import os

def catTool(process, runOnMC=True, useMiniAOD=True):
    if runOnMC:
        from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
        process.pileupWeight.pileupMC = pileupWeightMap[cat.pileupMCmap]
        process.pileupWeight.pileupRD = pileupWeightMap["%s"%cat.lumiJSON]
        process.pileupWeight.pileupUp = pileupWeightMap["%s_Up"%cat.lumiJSON]
        process.pileupWeight.pileupDn = pileupWeightMap["%s_Dn"%cat.lumiJSON]
        process.pileupWeightSilver = process.pileupWeight.clone()
        process.pileupWeightSilver.pileupRD = pileupWeightMap["%s"%cat.lumiJSONSilver]
        process.pileupWeightSilver.pileupUp = pileupWeightMap["%s_Up"%cat.lumiJSONSilver]
        process.pileupWeightSilver.pileupDn = pileupWeightMap["%s_Dn"%cat.lumiJSONSilver]
    else:
        from FWCore.PythonUtilities.LumiList import LumiList
        lumiJSON = os.environ["CMSSW_BASE"]+("/src/CATTools/CatProducer/data/LumiMask/%s.txt"%cat.lumiJSON)
        lumiJSONSilver = os.environ["CMSSW_BASE"]+("/src/CATTools/CatProducer/data/LumiMask/%s.txt"%cat.lumiJSONSilver)
        process.lumiMask = cms.EDFilter("LumiMaskFilter",
            LumiSections = LumiList(lumiJSON).getVLuminosityBlockRange())
        process.lumiMaskSilver = cms.EDFilter("LumiMaskFilter",
            LumiSections = LumiList(lumiJSONSilver).getVLuminosityBlockRange())
    
    useJECfile = True
    jecFile = cat.JetEnergyCorrection
    if runOnMC:
        jecFile = jecFile+"_MC"
    else:
        jecFile = jecFile+"_DATA"
    if useJECfile:
        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
            connect = cms.string('sqlite_fip:CATTools/CatProducer/data/JEC/%s.db'%jecFile),            
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_%s_AK4PF"%jecFile),
                    label= cms.untracked.string("AK4PF")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFchs"%jecFile),
                    label= cms.untracked.string("AK4PFchs")),
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFPuppi"%jecFile),
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
        ### updating puppi jet jec
        process.patJetPuppiCorrFactorsUpdated = process.patJetCorrFactorsUpdated.clone(
            src = process.catJetsPuppi.src,
            payload = cms.string('AK4PFPuppi'),
            levels = cms.vstring('L2Relative','L3Absolute'),
            useRho = cms.bool(False))
        
        process.patJetsPuppiUpdated = process.patJetsUpdated.clone(
            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetPuppiCorrFactorsUpdated")),
            jetSource = process.catJetsPuppi.src )
        ### updating pile Jet.
        process.load("RecoJets.JetProducers.PileupJetID_cfi")
        process.pileupJetIdUpdated = process.pileupJetId.clone(
          jets=cms.InputTag("slimmedJets"),
          inputIsCorrected=True,
          applyJec=True,
          vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
        )
        #process.patJetsUpdated.userData.userFloats.src +=['pileupJetIdUpdated:fullDiscriminant']

        process.catJets.src = cms.InputTag("patJetsUpdated")
        
        process.catJetsPuppi.src = cms.InputTag("patJetsPuppiUpdated")
        process.catJetsPuppi.setGenParticle = cms.bool(False)
        ## #######################################################################
        ## Setup JER
        ## JER needs random numbers
        process.RandomNumberGeneratorService.catJets = cms.PSet(
            engineName = cms.untracked.string('TRandom3'),
            initialSeed = cms.untracked.uint32(1),
        )
        process.RandomNumberGeneratorService.catJetsPuppi = cms.PSet(
            engineName = cms.untracked.string('TRandom3'),
            initialSeed = cms.untracked.uint32(1),
        )
        ## #######################################################################
        ## # MET corrections from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
        #from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
        #runMetCorAndUncFromMiniAOD( process, isData= not runOnMC, jecUncFile=cat.JECUncertaintyFile, jetColl= process.catJets.src)
        #process.catMETs.src = cms.InputTag("slimmedMETs","","CAT")
        #del process.slimmedMETs.caloMET
        ## redoing noHF met due to new correction
        #process.noHFCands = cms.EDFilter("CandPtrSelector",src=cms.InputTag("packedPFCandidates"),
        #                                 cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0"))
        #runMetCorAndUncFromMiniAOD(process,isData=not runOnMC,pfCandColl=cms.InputTag("noHFCands"),postfix="NoHF",
        #                           jecUncFile=cat.JECUncertaintyFile, jetColl= process.catJets.src)
        #process.catMETsNoHF = process.catMETs.clone(src = cms.InputTag("slimmedMETsNoHF","","CAT"))
        #del process.slimmedMETsNoHF.caloMET        
        #######################################################################
        ## for egamma pid https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_74X
        from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDPhotonIdProducer,setupAllVIDIdsInModule,setupVIDPhotonSelection            
        switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

        #Use IDs from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/RecoEgamma/PhotonIdentification/python/Identification/ 
        photon_ids = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                      'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']
        
        #add them to the VID producer
        for idmod in photon_ids:
            setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

        process.catPhotons.photonIDSources = cms.PSet( 
            cutBasedPhotonID_Spring15_25ns_V1_standalone_loose = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
            cutBasedPhotonID_Spring15_25ns_V1_standalone_medium = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
            cutBasedPhotonID_Spring15_25ns_V1_standalone_tight = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
            mvaPhoID_Spring15_25ns_nonTrig_V2_wp90 =  cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"),
            )
        
        ## from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection
        ## electron_ids = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
        ##                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
        ##                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
        ##                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']
        ## switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
        ## for idmod in electron_ids:
        ##     setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

        ## process.catElectrons.electronIDSources = cms.PSet(
        ##     cutBasedElectronID_Spring15_25ns_V1_standalone_loose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
        ##     cutBasedElectronID_Spring15_25ns_V1_standalone_medium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        ##     cutBasedElectronID_Spring15_25ns_V1_standalone_tight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
        ##     cutBasedElectronID_Spring15_25ns_V1_standalone_veto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
        ##     heepElectronID_HEEPV60 = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
        ##     mvaEleID_Spring15_25ns_nonTrig_V1_wp80 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
        ##     mvaEleID_Spring15_25ns_nonTrig_V1_wp90 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
        ##     mvaEleID_Spring15_25ns_Trig_V1_wp80 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
        ##     mvaEleID_Spring15_25ns_Trig_V1_wp90 = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
        ## )

        #######################################################################
        ## Energy smearing and scale correction
        ## https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
        process.RandomNumberGeneratorService.calibratedPatElectrons=cms.PSet(
            engineName = cms.untracked.string('TRandom3'),
            initialSeed = cms.untracked.uint32(1)
        )
        process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
        process.calibratedPatElectrons.isMC = runOnMC
        process.catElectrons.src = cms.InputTag("calibratedPatElectrons")    
    
        # photons not yet working...
        #process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
        
        #######################################################################    
        # adding pfMVAMet https://twiki.cern.ch/twiki/bin/viewauth/CMS/MVAMet#Spring15_samples_with_25ns_50ns
        # https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/RecoMET/METPUSubtraction/test/mvaMETOnMiniAOD_cfg.py
    ##     process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    ##     process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
    ##     process.ak4PFJets.doAreaFastjet = cms.bool(True)

    ## from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

    ## process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
    ## #process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons")
    ## process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
    ## process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")

    ## process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
    ## #process.puJetIdForPFMVAMEt.jets = cms.InputTag("ak4PFJets")
    ## process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
    ## process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")

    ## process.pfMVAMEt.inputFileNames = cms.PSet(
    ##     U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_25NS_July2015.root'),
    ##     DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_25NS_July2015.root'),
    ##     CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_25NS_July2015.root'),
    ##     CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_25NS_July2015.root')
    ## )
        
    ## process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
    ## process.patMETsPfMva = process.patMETs.clone(addGenMET = cms.bool(False), metSource  = cms.InputTag("pfMVAMEt"))
    ## process.catMETsPfMva = process.catMETs.clone(src = cms.InputTag("patMETsPfMva"))
    ## process.catMETsPfMva.setUnclusteredEn = cms.bool(False)
    ## process.catMETsPfMva.setJetMETSyst = cms.bool(False)
    
