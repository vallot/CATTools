dataset_loc=$CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset
cfg=$CMSSW_BASE/src/CATTools/CatAnalyzer/prod/ttbbLepJetsAnalyzer_cfg.py
save_loc=/store/user/san/ntuple/Run2018/V10_2

# 2017 Run2 data
### EG
create-batch --jobName DataSingleEGA --fileList $dataset_loc/dataset_EGamma_Run2018A.txt --cfg $cfg --transferDest $save_loc/DataSingleEGA --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGB --fileList $dataset_loc/dataset_EGamma_Run2018B.txt --cfg $cfg --transferDest $save_loc/DataSingleEGB --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGC --fileList $dataset_loc/dataset_EGamma_Run2018C.txt --cfg $cfg --transferDest $save_loc/DataSingleEGC --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGD --fileList $dataset_loc/dataset_EGamma_Run2018D.txt --cfg $cfg --transferDest $save_loc/DataSingleEGD --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
### Mu
create-batch --jobName DataSingleMuA --fileList $dataset_loc/dataset_SingleMuon_Run2018A.txt --cfg $cfg --transferDest $save_loc/DataSingleMuB --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuB --fileList $dataset_loc/dataset_SingleMuon_Run2018B.txt --cfg $cfg --transferDest $save_loc/DataSingleMuB --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuC --fileList $dataset_loc/dataset_SingleMuon_Run2018C.txt --cfg $cfg --transferDest $save_loc/DataSingleMuC --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuD --fileList $dataset_loc/dataset_SingleMuon_Run2018D.txt --cfg $cfg --transferDest $save_loc/DataSingleMuD --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

# Signal
create-batch --jobName TTLJ_PowhegPythia_ttbb --fileList $dataset_loc/dataset_TTLJ_powheg.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_ttbb --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'
create-batch --jobName TTLJ_PowhegPythia_ttbj --fileList $dataset_loc/dataset_TTLJ_powheg.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_ttbj --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'
create-batch --jobName TTLJ_PowhegPythia_ttcc --fileList $dataset_loc/dataset_TTLJ_powheg.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_ttcc --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'
create-batch --jobName TTLJ_PowhegPythia_ttLF --fileList $dataset_loc/dataset_TTLJ_powheg.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_ttLF --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'
create-batch --jobName TTLJ_PowhegPythia_ttother --fileList $dataset_loc/dataset_TTLJ_powheg.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_ttother --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=5'

# Background
### ttbar background
create-batch --jobName TTLJ_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTLJ_powheg.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythiaBkg --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTLL_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTTo2L2Nu.txt --cfg $cfg --transferDest $save_loc/TTLL_PowhegPythiaBkg --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTJJ_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTToHadronic.txt --cfg $cfg --transferDest $save_loc/TTJJ_PowhegPythiaBkg --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
### single top
create-batch --jobName SingleTop_s_aMCatNLOPythia --fileList $dataset_loc/dataset_SingleTop_s.txt --cfg $cfg --transferDest $save_loc/SingleTop_s_aMCatNLOPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName SingleTop_t_PowhegPythia --fileList $dataset_loc/dataset_SingleTop_t.txt --cfg $cfg --transferDest $save_loc/SingleTop_t_PowhegPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName SingleTop_tW_PowhegPythia --fileList $dataset_loc/dataset_SingleTop_tW.txt --cfg $cfg --transferDest $save_loc/SingleTop_tW_PowhegPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName SingleTbar_t_PowhegPythia --fileList $dataset_loc/dataset_SingleTbar_t.txt --cfg $cfg --transferDest $save_loc/SingleTbar_t_PowhegPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName SingleTbar_tW_PowhegPythia --fileList $dataset_loc/dataset_SingleTbar_tW.txt --cfg $cfg --transferDest $save_loc/SingleTbar_tW_PowhegPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
### ttbar + boson
create-batch --jobName ttHToNonbb_PowhegPythia --fileList $dataset_loc/dataset_ttHToNonbb.txt --cfg $cfg --transferDest $save_loc/ttHToNonbb_PowhegPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName ttHTobb_PowhegPythia --fileList $dataset_loc/dataset_ttHTobb.txt --cfg $cfg --transferDest $save_loc/ttHTobb_PowhegPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName ttWToLNu_aMCatNLOMadspinPythia --fileList $dataset_loc/dataset_TTWJetsToLNu_PSweight.txt --cfg $cfg --transferDest $save_loc/ttWToLNu_aMCatNLOMadspinPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName ttWToQQ_aMCatNLOMadspinPythia --fileList $dataset_loc/dataset_TTWJetsToQQ.txt --cfg $cfg --transferDest $save_loc/ttWToQQ_aMCatNLOMadspinPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName ttZToLLNuNu_aMCatNLOMadspinPythia --fileList $dataset_loc/dataset_TTZToLLNuNu.txt --cfg $cfg --transferDest $save_loc/ttZToLLNuNu_aMCatNLOMadspinPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName ttZToQQ_aMCatNLOMadspinPythia --fileList $dataset_loc/dataset_TTZToQQ.txt --cfg $cfg --transferDest $save_loc/ttZToQQ_aMCatNLOMadspinPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
### diboson
create-batch --jobName WW_Pythia  --fileList $dataset_loc/dataset_WW.txt --cfg $cfg --transferDest $save_loc/WW_Pythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName WZ_Pythia  --fileList $dataset_loc/dataset_WZ.txt --cfg $cfg --transferDest $save_loc/WZ_Pythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName ZZ_Pythia  --fileList $dataset_loc/dataset_ZZ.txt --cfg $cfg --transferDest $save_loc/ZZ_Pythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
### W + Jets
create-batch --jobName WJetsToLNu_MadgraphPythia  --fileList $dataset_loc/dataset_WJetsToLNu.txt --cfg $cfg --transferDest $save_loc/WJetsToLNu_MadgraphPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName WJetsToLNu_MadgraphPythia_part2  --fileList $dataset_loc/dataset_WJetsToLNu_part2.txt --cfg $cfg --transferDest $save_loc/WJetsToLNu_MadgraphPythia_part2 --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName W1JetsToLNu_MadgraphPythia  --fileList $dataset_loc/dataset_W1JetsToLNu.txt --cfg $cfg --transferDest $save_loc/W1JetsToLNu_MadgraphPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName W2JetsToLNu_MadgraphPythia  --fileList $dataset_loc/dataset_W2JetsToLNu.txt --cfg $cfg --transferDest $save_loc/W2JetsToLNu_MadgraphPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName W3JetsToLNu_MadgraphPythia  --fileList $dataset_loc/dataset_W3JetsToLNu.txt --cfg $cfg --transferDest $save_loc/W3JetsToLNu_MadgraphPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName W4JetsToLNu_MadgraphPythia  --fileList $dataset_loc/dataset_W4JetsToLNu.txt --cfg $cfg --transferDest $save_loc/W4JetsToLNu_MadgraphPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
### Z + Jets
create-batch --jobName ZJets_M10to50_MadgraphPythia  --fileList $dataset_loc/dataset_DYJets_10to50.txt --cfg $cfg --transferDest $save_loc/ZJets_M10to50_MadgraphPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName ZJets_M50_aMCatNLOPythia  --fileList $dataset_loc/dataset_DYJets.txt --cfg $cfg --transferDest $save_loc/ZJets_M50_aMCatNLOPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
