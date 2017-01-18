import FWCore.ParameterSet.Config as cms

## Energy smearing and scale correction
## https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression

def enableElectronRegression(process):
    from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
    process = regressionWeights(process)

    process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

    return process

