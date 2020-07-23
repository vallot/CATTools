import FWCore.ParameterSet.Config as cms
import catDefinitions_cfi as cat
import os
print os.environ['CMSSW_BASE']

def catTool(process, runOnMC=True, useMiniAOD=True):
    if runOnMC:
        from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
        process.pileupWeight.pileupMC = pileupWeightMap[cat.pileupMCmap]
        process.pileupWeight.pileupRD = pileupWeightMap["%s"%cat.lumiJSON]
        process.pileupWeight.pileupUp = pileupWeightMap["%s_Up"%cat.lumiJSON]
        process.pileupWeight.pileupDn = pileupWeightMap["%s_Dn"%cat.lumiJSON]
    else:
        from FWCore.PythonUtilities.LumiList import LumiList
        process.lumiMask = cms.EDFilter("LumiMaskFilter",
            LumiSections = LumiList('%s/src/CATTools/CatProducer/data/LumiMask/%s.txt'%(os.environ['CMSSW_BASE'], cat.lumiJSON)).getVLuminosityBlockRange())

        #process.load("CATTools.CatProducer.eventCleaning.badECALSlewRateMitigationFilter2016_cfi")

    useJECfile = True
    jecFiles = cat.JetEnergyCorrection
    if runOnMC:
        jecFile = jecFiles[1]
    else:
        jecFile = jecFiles[0]
    if useJECfile:
        from CondCore.CondDB.CondDB_cfi import CondDB
        if hasattr(CondDB, 'connect'): delattr(CondDB, 'connect')
        process.jec = cms.ESSource("PoolDBESSource",CondDB,
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

        ## applying new jec on the fly
#        process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
#        process.updatedPatJetCorrFactors.levels = cms.vstring('L1FastJet','L2Relative','L3Absolute','L2L3Residual')
#        process.patJetCorrFactors.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
#        process.catJets.src = cms.InputTag("updatedPatJets")

        from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
        updateJetCollection(
           process,
           jetSource = cms.InputTag('slimmedJets'),
           pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
           svSource = cms.InputTag('slimmedSecondaryVertices'),
           jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
           btagDiscriminators = [
              'pfDeepFlavourJetTags:probb',
              'pfDeepFlavourJetTags:probbb',
              'pfDeepFlavourJetTags:problepb',
              'pfDeepFlavourJetTags:probc',
              'pfDeepFlavourJetTags:probuds',
              'pfDeepFlavourJetTags:probg'
              ],
           postfix='NewDFTraining'
        )
        process.p += process.patJetCorrFactorsNewDFTraining
        process.p += process.updatedPatJetsNewDFTraining

        updateJetCollection(
           process,
           jetSource = cms.InputTag('slimmedJets'),
           pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
           svSource = cms.InputTag('slimmedSecondaryVertices'),
           labelName = 'UpdatedJEC',
           jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
        )
        process.p += process.patJetCorrFactorsUpdatedJEC
        process.p += process.updatedPatJetsUpdatedJEC
        process.catJets.src = cms.InputTag("updatedPatJetsUpdatedJEC","","CAT")

#    if useMiniAOD: ## corrections when using miniAOD #This is stored in miniAOD as flag, in 2017
#        from CATTools.CatProducer.patTools.metFilters_cff import enableAdditionalMETFilters
#        process = enableAdditionalMETFilters(process, runOnMC)

        #######################################################################
        # puppi https://twiki.cern.ch/twiki/bin/view/CMS/PUPPI
        # using default
        #######################################################################
        ### updating puppi jet jec
#        process.patJetPuppiCorrFactorsUpdated = process.updatedPatJetCorrFactors.clone(
#            src = process.catJetsPuppi.src,
#            payload = cms.string('AK4PFPuppi'),
#            levels = cms.vstring('L2Relative','L3Absolute'),
#            useRho = cms.bool(False))
        
#        process.patJetsPuppiUpdated = process.updatedPatJets.clone(
#            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetPuppiCorrFactorsUpdated")),
#            jetSource = process.catJetsPuppi.src )
        ### updating pile Jet.
#        process.load("RecoJets.JetProducers.PileupJetID_cfi")
#        process.pileupJetIdUpdated = process.pileupJetId.clone(
#          jets=cms.InputTag("slimmedJets"),
#          inputIsCorrected=True,
#          applyJec=True,
#          vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
#        )
        #process.patJetsUpdated.userData.userFloats.src +=['pileupJetIdUpdated:fullDiscriminant']

#        process.catJetsPuppi.src = cms.InputTag("patJetsPuppiUpdated")
#        process.catJetsPuppi.setGenParticle = cms.bool(False)
        ## #######################################################################
        ## qg-likelihood
        # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
#        from CATTools.CatProducer.patTools.jetQGLikelihood_cff import enableQGLikelihood
#        process = enableQGLikelihood(process, qgDatabaseVersion="v2b", runOnMC=runOnMC, useMiniAOD=useMiniAOD)

        ## DeepFlavour
        #from CATTools.CatProducer.patTools.jetDeepFlavour_cff import enableDeepFlavour
        #process = enableDeepFlavour(process)

        ## #######################################################################
        ## # MET corrections from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
#       from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#       runMetCorAndUncFromMiniAOD(process, isData= not runOnMC, electronColl=cms.InputTag('calibratedPatElectrons'))
#       process.catMETs.src = cms.InputTag("slimmedMETs","","CAT")

        #from CATTools.CatProducer.patTools.metMuonRecoMitigation2016_cff import enableMETMuonRecoMitigation2016
        #process = enableMETMuonRecoMitigation2016(process, runOnMC) ## MET input object is overridden in the modifier function

    if useMiniAOD:
      # Instructions for 9_4_X, X >=9 for 2017 data with EE noise mitigation
      from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

      runMetCorAndUncFromMiniAOD (
              process,
              isData = not runOnMC,
              fixEE2017 = True,
              fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139},
              postfix = "ModifiedMET"
      )

      process.p += process.fullPatMetSequenceModifiedMET
      process.catMETs.src = cms.InputTag("slimmedMETsModifiedMET","","CAT")

      # https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
      from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
      process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
          DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
          UseJetEMPt = cms.bool(False),
          PrefiringRateSystematicUncty = cms.double(0.2),
          SkipWarnings = False)

      process.p += process.prefiringweight

      # Run ecalBadCalibFilter: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
      process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
      baddetEcallist = cms.vuint32([872439604,872422825,872420274,872423218,
           872423215,872416066,872435036,872439336,872420273,872436907,872420147,872439731,
           872436657,872420397,872439732,872439339,872439603,872422436,872439861,872437051,
           872437052,872420649,872422436,872421950,872437185,872422564,872421566,872421695,
           872421955,872421567,872437184,872421951,872421694,872437056,872437057,872437313])
      process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
          "EcalBadCalibFilter",
          EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
          ecalMinEt        = cms.double(50.),
          baddetEcal    = baddetEcallist,
          taggingMode = cms.bool(True),
          debug = cms.bool(False))
      process.p += process.ecalBadCalibReducedMINIAODFilter

      # Egamma post reco tools: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#Running_on_2017_MiniAOD_V2
      from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
      setupEgammaPostRecoSeq(process,
                             runVID=False, #if you want the Fall17V2 IDs, set this to True or remove (default is True)
                             era='2017-Nov17ReReco')  #era is new to select between 2016 / 2017,  it defaults to 2017
      process.p += process.egammaPostRecoSeq

def addEgmID(process, runOnMC):
        #######################################################################
        ## Electron regression
        #from CATTools.CatProducer.patTools.egmRegression_cff import enableElectronRegression
        #process = enableElectronRegression(process)

        ## Energy/Photon smearing and scale correction
        #from CATTools.CatProducer.patTools.egmSmearing_cff import enableElectronSmearing, enablePhotonSmearing
        #process = enableElectronSmearing(process, runOnMC)
        #process = enablePhotonSmearing(process, runOnMC)
        
        ## Electron/Photon VID
        from CATTools.CatProducer.patTools.egmVersionedID_cff import enableElectronVID, enablePhotonVID
        process = enableElectronVID(process)
        process.p += process.egmGsfElectronIDSequence
        #process = enablePhotonVID(process)

        ## Electron ID without isolation cuts
        from CATTools.CatProducer.patTools.egmNoIsoID_cff import enableElectronNoIsoID
        process = enableElectronNoIsoID(process)
