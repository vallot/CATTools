import FWCore.ParameterSet.Config as cms
# Jet quark-gluon likelihood
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

from CondCore.DBCommon.CondDBSetup_cfi import *
def enableQGLikelihood(process, qgDatabaseVersion='v2b', runOnMC=True, useMiniAOD=True):

    process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(),
        #connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000')
        connect = cms.string('sqlite_fip:CATTools/CatProducer/data/QGL_%s.db' % qgDatabaseVersion),
    )

    process.load('RecoJets.JetProducers.QGTagger_cfi')
    for type in ['AK4PFchs']:#, 'AK4PF',]:
        process.QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
            label  = cms.untracked.string('QGL_'+type)
        )))
    process.es_prefer_QGPoolDBESSource = cms.ESPrefer("PoolDBESSource","QGPoolDBESSource")

    process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')
    #if useMiniAOD:
    #    process.QGTagger.srcJets = cms.InputTag("slimmedJets")
    #else:
    #    process.QGTagger.srcJets = cms.InputTag("updatedPatJets")
    #    #process.QGTagger.srcJets    = cms.InputTag("selectedPatJetsAK4PFCHS")

    process.catJets.qgLikelihood = 'QGTagger:qgLikelihood'
    process.QGTagger.srcJets = process.catJets.src

    return process
