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
            )
        )
        process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
        print "JEC based on", process.jec.connect

        ## applying new jec on the fly
        from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
        updateJetCollection(
           process,
           jetSource = cms.InputTag('slimmedJets'),
           labelName = 'UpdatedJEC',
           jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
           #pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
           #svSource = cms.InputTag('slimmedSecondaryVertices'),
           #btagDiscriminators = [ #included in MiniAODv3
           #   'pfDeepFlavourJetTags:probb',
           #   'pfDeepFlavourJetTags:probbb',
           #   'pfDeepFlavourJetTags:problepb',
           #   'pfDeepFlavourJetTags:probc',
           #   'pfDeepFlavourJetTags:probuds',
           #   'pfDeepFlavourJetTags:probg'
           #   ],
        )
        process.p += process.patJetCorrFactorsUpdatedJEC
        process.p += process.updatedPatJetsUpdatedJEC
        process.catJets.src = cms.InputTag("updatedPatJetsUpdatedJEC","","CAT")

        #######################################################################
        # puppi https://twiki.cern.ch/twiki/bin/view/CMS/PUPPI
        # using default
        #######################################################################
        ### updating pile Jet.
#        process.load("RecoJets.JetProducers.PileupJetID_cfi")
#        process.pileupJetIdUpdated = process.pileupJetId.clone(
#          jets=cms.InputTag("slimmedJets"),
#          inputIsCorrected=True,
#          applyJec=True,
#          vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
#        )
        #process.patJetsUpdated.userData.userFloats.src +=['pileupJetIdUpdated:fullDiscriminant']


    ## #######################################################################
    ## # MET corrections from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
    if useMiniAOD:
      from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

      runMetCorAndUncFromMiniAOD (
              process,
              isData = not runOnMC,
              postfix = "ModifiedMET"
      )

      process.p += process.fullPatMetSequenceModifiedMET
      process.catMETs.src = cms.InputTag("slimmedMETsModifiedMET","","CAT")

      # https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
      from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
      process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
          DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
          UseJetEMPt = cms.bool(False),
          PrefiringRateSystematicUncty = cms.double(0.2),
          SkipWarnings = False)
      process.p += process.prefiringweight

      # Egamma post reco tools: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#Running_on_2016_MiniAOD_V2
      from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
      setupEgammaPostRecoSeq(process,
                             runVID=False, #run in old way, in egmVersionedID_cff.py
                             runEnergyCorrections=False,
                             era='2016-Legacy')
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
