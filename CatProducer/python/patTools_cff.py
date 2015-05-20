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

        ## tem - due to ak4GenJetsNoNu not in AOD for now
        process.load('CommonTools.ParticleFlow.genForPF2PAT_cff')
