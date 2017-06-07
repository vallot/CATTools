import FWCore.ParameterSet.Config as cms

commonTestCATTuples = {
    "sig":cms.untracked.vstring("/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_103557/0000/catTuple_1.root"),
    "bkg":cms.untracked.vstring("/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_102707/0000/catTuple_1.root"),
    "data":cms.untracked.vstring("/store/group/CAT/DoubleEG/v8-0-6_Run2016H-03Feb2017_ver3-v1/170303_095044/0000/catTuple_1.root"),
}

commonTestCATTuples804 = {
    "sig":cms.untracked.vstring("/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-4_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170113_045423/0000/catTuple_1.root"),
    "bkg":cms.untracked.vstring("/store/group/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/v8-0-4_RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170116_095242/0000/catTuple_1.root"),
    "data":cms.untracked.vstring("/store/group/CAT/DoubleEG/v8-0-4_Run2016H-PromptReco-v3/170113_133053/0000/catTuple_1.root"),
}

commonTestMiniAODs = {
    "sig":cms.untracked.vstring("root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root",),
    "bkg":cms.untracked.vstring("root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0AF0207B-EFBE-E611-B4BE-0CC47A7FC858.root",),
    "data":cms.untracked.vstring("root://cmsxrootd.fnal.gov///store/data/Run2016H/MuonEG/MINIAOD/03Feb2017_ver2-v1/100000/044366C7-4AEE-E611-8CF7-0025905B856E.root",),
}
