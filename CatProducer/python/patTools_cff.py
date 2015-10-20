import FWCore.ParameterSet.Config as cms

def patTool(process, runOnMC=True, useMiniAOD = True):
    if not useMiniAOD:#running pretty much the default miniAOD steps!
        # doing all filters
        process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
        from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import allMetFilterPaths
        for filt in allMetFilterPaths:
            process.schedule.append(getattr(process,'Flag_'+filt))

        process.load('Configuration.StandardSequences.PAT_cff')
        
        from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC,miniAOD_customizeAllData
        if runOnMC:
            miniAOD_customizeAllMC(process)
        else :
            miniAOD_customizeAllData(process)
        
        process.patJetsPuppi.embedGenPartonMatch = cms.bool(False)
        process.patJetCorrFactorsPuppi.useRho = cms.bool(False)
    process.schedule.append(process.p)
