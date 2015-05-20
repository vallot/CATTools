import FWCore.ParameterSet.Config as cms

def patTool(process, runOnMC=True, useMiniAOD = True):
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")

    if not useMiniAOD:
        process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
        process.load('Configuration.StandardSequences.PAT_cff')
        from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC,miniAOD_customizeAllData
        if runOnMC:
            miniAOD_customizeAllMC(process)
        else :
            miniAOD_customizeAllData(process)

        #for muon isolation
        #process.patMuons.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03","muPFIsoValueNeutral03","muPFIsoValueGamma03","muPFIsoValuePU03","muPFIsoValueChargedAll03")
        #for electron isolation
        process.patElectrons.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFId","elPFIsoValueNeutral03PFId","elPFIsoValueGamma03PFId","elPFIsoValuePU03PFId","elPFIsoValueChargedAll03PFId")
        process.patElectrons.isolationValuesNoPFId.user = cms.VInputTag("elPFIsoValueCharged03NoPFId","elPFIsoValueNeutral03NoPFId","elPFIsoValueGamma03NoPFId","elPFIsoValuePU03NoPFId","elPFIsoValueChargedAll03NoPFId")

        
