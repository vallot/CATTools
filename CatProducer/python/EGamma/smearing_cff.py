import FWCore.ParameterSet.Config as cms

## Energy smearing and scale correction
## https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer

def enableElectronSmearing(process, runOnMC=True):
    process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

    process.RandomNumberGeneratorService.calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(81),
        engineName = cms.untracked.string('TRandom3')
    )

    process.calibratedPatElectrons.isMC = runOnMC
    process.catElectrons.src = "calibratedPatElectrons"

    return process

def enablePhotonSmearing(process, runOnMC=True):
    process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')

    process.RandomNumberGeneratorService.calibratedPatPhotons = cms.PSet(
        initialSeed = cms.untracked.uint32(81),
        engineName = cms.untracked.string('TRandom3')
    )

    process.calibratedPatPhotons.isMC = runOnMC
    process.catPhotons.src = "calibratedPatPhotons"

    return process

