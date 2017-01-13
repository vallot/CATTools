import FWCore.ParameterSet.Config as cms
# Jet quark-gluon likelihood
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

from CondCore.DBCommon.CondDBSetup_cfi import *
def enableDeepFlavour(process):
    process.catJets.flavTagLabels.extend([
        cms.InputTag("deepFlavourJetTags:probudsg"),
        cms.InputTag("deepFlavourJetTags:probb"),
        cms.InputTag("deepFlavourJetTags:probc"),
        cms.InputTag("deepFlavourJetTags:probbb"),
        cms.InputTag("deepFlavourJetTags:probcc"),
    ])

    return process
