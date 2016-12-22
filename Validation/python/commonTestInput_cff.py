import FWCore.ParameterSet.Config as cms

commonTestCATTuples = {
    "sig":cms.untracked.vstring("/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v8-0-3_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/161205_142558/0000/catTuple_1.root",),
    "bkg":cms.untracked.vstring("/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-3_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/161205_142228/0000/catTuple_1.root",),
    "data":cms.untracked.vstring("/store/group/CAT/DoubleEG/v8-0-3_Run2016G-23Sep2016-v1/161204_005359/0000/catTuple_1.root"),
}

commonTestMiniAODs = {
    "sig":cms.untracked.vstring("root://cmsxrootd.fnal.gov//store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/004A0552-3929-E611-BD44-0025905A48F0.root",),
    "bkg":cms.untracked.vstring("root://cmsxrootd.fnal.gov//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root",),
    "data":cms.untracked.vstring("root://cmsxrootd.fnal.gov///store/data/Run2016C/SingleMuon/MINIAOD/23Sep2016-v1/80000/82AF08FC-3B87-E611-A209-FA163E3F4268.root",),
}
