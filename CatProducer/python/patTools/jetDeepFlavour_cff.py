import FWCore.ParameterSet.Config as cms
# https://twiki.cern.ch/twiki/bin/view/CMS/DeepFlavour
from CondCore.DBCommon.CondDBSetup_cfi import *
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

def enableDeepFlavour(process):
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = ["pfDeepCSVJetTags:probudsg",
                              "pfDeepCSVJetTags:probb",
                              "pfDeepCSVJetTags:probc",
                              "pfDeepCSVJetTags:probbb",
                              "pfDeepCSVJetTags:probcc",
                              "pfDeepCMVAJetTags:probudsg",
                              "pfDeepCMVAJetTags:probb",
                              "pfDeepCMVAJetTags:probc",
                              "pfDeepCMVAJetTags:probbb",
                              "pfDeepCMVAJetTags:probcc",
                             ]
    )

    process.catJets.src = "updatedPatJetsTransientCorrected"

    return process
