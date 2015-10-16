import FWCore.ParameterSet.Config as cms

def patTool(process, runOnMC=True, useMiniAOD = True):
    if not useMiniAOD:
        # doing all filters
        process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
        from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import allMetFilterPaths
        for filt in allMetFilterPaths:
            process.schedule.append(getattr(process,'Flag_'+filt))

        process.load('Configuration.StandardSequences.PAT_cff')
        #process.p += process.selectedPatTrigger + process.slimmedElectrons + process.slimmedMuons + process.slimmedJets + process.slimmedMETs
        process.schedule.append(process.p)
       
        from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC,miniAOD_customizeAllData
        if runOnMC:
            miniAOD_customizeAllMC(process)
        else :
            miniAOD_customizeAllData(process)

        # because of unscheduled, for miniAOD to run first
        #process.p += process.slimmedJets + process.slimmedMuons + process.slimmedElectrons
