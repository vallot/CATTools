import FWCore.ParameterSet.Config as cms

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Recipe_details_80X_94X
prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
    ThePhotons = cms.InputTag("slimmedPhotons"),
    TheJets = cms.InputTag("slimmedJets"),
    L1Maps = cms.string("../../../L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root"), # update this line with the location of this file
    DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
    PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
)
