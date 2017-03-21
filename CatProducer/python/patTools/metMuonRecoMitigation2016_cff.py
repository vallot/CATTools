import FWCore.ParameterSet.Config as cms

def enableMETMuonRecoMitigation2016(process, runOnMC):
    ## using updated recipes from https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes
    # Now you are creating the bad muon corrected MET
    from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation
    muonRecoMitigation(
        process = process,
        pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
        runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
        selection="", #You can use a custom selection for your bad muons. Leave empty if you would like to use the bad muon recipe definition.
        muonCollection="", #The muon collection name where your custom selection will be applied to. Leave empty if you would like to use the bad muon recipe definition.
        cleanCollName="cleanMuonsPFCandidates", #output pf candidate collection ame
        cleaningScheme="computeAllApplyClone", #Options are: "all", "computeAllApplyBad","computeAllApplyClone". Decides which (or both) bad muon collections to be used for MET cleaning coming from the bad muon recipe.
        postfix="" #Use if you would like to add a post fix to your muon / pf collections
    )

    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process,
        isData= not runOnMC,
        pfCandColl="cleanMuonsPFCandidates",
        recoMetFromPFCs=True,
        postfix="MuClean"
    )
     
    process.mucorMET = cms.Sequence(                     
        process.badGlobalMuonTaggerMAOD *
        process.cloneGlobalMuonTaggerMAOD *
        #process.badMuons * # If you are using cleaning mode "all", uncomment this line
        process.cleanMuonsPFCandidates *
        process.fullPatMetSequenceMuClean
    )

    process.catMETs.src = cms.InputTag("slimmedMETsMuClean","","CAT")                                                                                                          

    return process
