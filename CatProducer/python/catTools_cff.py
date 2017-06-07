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

        process.load("CATTools.CatProducer.eventCleaning.badECALSlewRateMitigationFilter2016_cfi")
    
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
    
    if useMiniAOD: ## corrections when using miniAOD
        from CATTools.CatProducer.patTools.metFilters_cff import enableAdditionalMETFilters
        process = enableAdditionalMETFilters(process, runOnMC)

        #######################################################################
        # puppi https://twiki.cern.ch/twiki/bin/view/CMS/PUPPI
        # using default
        #######################################################################
        ## applying new jec on the fly
        process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
        if not runOnMC:
            process.updatedPatJetCorrFactors.levels = cms.vstring('L1FastJet','L2Relative','L3Absolute','L2L3Residual')
            
        process.patJetCorrFactors.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
        process.catJets.src = cms.InputTag("updatedPatJets")
        ### updating puppi jet jec
        process.patJetPuppiCorrFactorsUpdated = process.updatedPatJetCorrFactors.clone(
            src = process.catJetsPuppi.src,
            payload = cms.string('AK4PFPuppi'),
            levels = cms.vstring('L2Relative','L3Absolute'),
            useRho = cms.bool(False))
        
        process.patJetsPuppiUpdated = process.updatedPatJets.clone(
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

        process.catJetsPuppi.src = cms.InputTag("patJetsPuppiUpdated")
        process.catJetsPuppi.setGenParticle = cms.bool(False)
        ## #######################################################################
        ## Setup JER
        ## JER needs random numbers
        process.RandomNumberGeneratorService.catJets = cms.PSet(
            engineName = cms.untracked.string('TRandom3'),
            initialSeed = cms.untracked.uint32(1),
        )
        process.RandomNumberGeneratorService.catFatJets = cms.PSet(
            engineName = cms.untracked.string('TRandom3'),
            initialSeed = cms.untracked.uint32(1),
        )
        process.RandomNumberGeneratorService.catJetsPuppi = cms.PSet(
            engineName = cms.untracked.string('TRandom3'),
            initialSeed = cms.untracked.uint32(1),
        )

        ## qg-likelihood
        # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
        from CATTools.CatProducer.patTools.jetQGLikelihood_cff import enableQGLikelihood
        process = enableQGLikelihood(process, qgDatabaseVersion="v2b", runOnMC=runOnMC, useMiniAOD=useMiniAOD)

        ## DeepFlavour
        from CATTools.CatProducer.patTools.jetDeepFlavour_cff import enableDeepFlavour
        process = enableDeepFlavour(process)

        ## #######################################################################
        ## # MET corrections from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
        runMetCorAndUncFromMiniAOD(process, isData= not runOnMC, electronColl=cms.InputTag('calibratedPatElectrons'))
        process.catMETs.src = cms.InputTag("slimmedMETs","","CAT")

        from CATTools.CatProducer.patTools.metMuonRecoMitigation2016_cff import enableMETMuonRecoMitigation2016
        process = enableMETMuonRecoMitigation2016(process, runOnMC) ## MET input object is overridden in the modifier function

        #######################################################################
        ## Electron regression
        from CATTools.CatProducer.patTools.egmRegression_cff import enableElectronRegression
        process = enableElectronRegression(process)

        ## Energy/Photon smearing and scale correction
        from CATTools.CatProducer.patTools.egmSmearing_cff import enableElectronSmearing, enablePhotonSmearing
        process = enableElectronSmearing(process, runOnMC)
        process = enablePhotonSmearing(process, runOnMC)
        
        ## Electron/Photon VID
        from CATTools.CatProducer.patTools.egmVersionedID_cff import enableElectronVID, enablePhotonVID
        process = enableElectronVID(process)
        process = enablePhotonVID(process)

        ## Electron ID without isolation cuts
        from CATTools.CatProducer.patTools.egmNoIsoID_cff import enableElectronNoIsoID
        process = enableElectronNoIsoID(process)
       
