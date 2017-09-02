import FWCore.ParameterSet.Config as cms
# https://twiki.cern.ch/twiki/bin/view/CMS/DeepFlavour

from CondCore.DBCommon.CondDBSetup_cfi import *
def enableDeepFlavour(process):
    process.load("RecoBTag.Combined.deepFlavour_cff")

    process.catJets.flavTagLabels.extend([
        cms.InputTag("pfDeepFlavourJetTags:probudsg"),
        cms.InputTag("pfDeepFlavourJetTags:probb"),
        cms.InputTag("pfDeepFlavourJetTags:probc"),
        cms.InputTag("pfDeepFlavourJetTags:probbb"),
        cms.InputTag("pfDeepFlavourJetTags:probcc"),
        cms.InputTag("pfDeepFlavourCMVAJetTags:probudsg"),
        cms.InputTag("pfDeepFlavourCMVAJetTags:probb"),
        cms.InputTag("pfDeepFlavourCMVAJetTags:probc"),
        cms.InputTag("pfDeepFlavourCMVAJetTags:probbb"),
        cms.InputTag("pfDeepFlavourCMVAJetTags:probcc"),
    ])

    return process
