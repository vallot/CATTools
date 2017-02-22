import FWCore.ParameterSet.Config as cms

## Energy smearing and scale correction
## https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression

def enableElectronRegression(process):
    from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
    process = regressionWeights(process)

    from EgammaAnalysis.ElectronTools.regressionWeights_local_cfi import GBRDWrapperRcd
    GBRDWrapperRcd.connect = "sqlite_fip:EgammaAnalysis/ElectronTools/data/ged_regression_20170114.db"

    process.regressions           = GBRDWrapperRcd
    process.es_prefer_regressions = cms.ESPrefer('PoolDBESSource','regressions')

    process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

    process.catElectrons.unsmaredElectrons = "slimmedElectrons::"+process.process

    return process

