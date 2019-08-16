dataset_loc=$CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset
cfg=$CMSSW_BASE/src/CATTools/CatAnalyzer/prod/ttbbLepJetsAnalyzer_cfg.py
save_loc=/store/user/san/ntuple/Run2017/V9_6

# 2017 Run2 data
### EG
create-batch --jobName DataSingleEGB --fileList $dataset_loc/dataset_SingleElectron_Run2017B.txt --cfg $cfg --transferDest $save_loc/DataSingleEGB --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGC --fileList $dataset_loc/dataset_SingleElectron_Run2017C.txt --cfg $cfg --transferDest $save_loc/DataSingleEGC --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGD --fileList $dataset_loc/dataset_SingleElectron_Run2017D.txt --cfg $cfg --transferDest $save_loc/DataSingleEGD --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGE --fileList $dataset_loc/dataset_SingleElectron_Run2017E.txt --cfg $cfg --transferDest $save_loc/DataSingleEGE --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleEGF --fileList $dataset_loc/dataset_SingleElectron_Run2017F.txt --cfg $cfg --transferDest $save_loc/DataSingleEGF --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
### Mu
create-batch --jobName DataSingleMuB --fileList $dataset_loc/dataset_SingleMuon_Run2017B.txt --cfg $cfg --transferDest $save_loc/DataSingleMuB --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuC --fileList $dataset_loc/dataset_SingleMuon_Run2017C.txt --cfg $cfg --transferDest $save_loc/DataSingleMuC --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuD --fileList $dataset_loc/dataset_SingleMuon_Run2017D.txt --cfg $cfg --transferDest $save_loc/DataSingleMuD --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuE --fileList $dataset_loc/dataset_SingleMuon_Run2017E.txt --cfg $cfg --transferDest $save_loc/DataSingleMuE --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName DataSingleMuF --fileList $dataset_loc/dataset_SingleMuon_Run2017F.txt --cfg $cfg --transferDest $save_loc/DataSingleMuF --maxFiles 50 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

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
create-batch --jobName ttWToLNu_aMCatNLOMadspinPythia --fileList $dataset_loc/dataset_TTWJetsToLNu.txt --cfg $cfg --transferDest $save_loc/ttWToLNu_aMCatNLOMadspinPythia --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
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

#QCD
### muon
create-batch --jobName QCD_Pt-1000toInf_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-1000toInf_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-1000toInf_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-800to1000_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-800to1000_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-800to1000_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-600to800_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-600to800_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-600to800_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-470to600_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-470to600_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-470to600_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-300to470_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-300to470_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-300to470_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-170to300_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-170to300_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-170to300_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-120to170_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-120to170_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-120to170_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-80to120_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-80to120_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-80to120_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-50to80_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-50to80_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-50to80_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-30to50_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-30ot050_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-30to50_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-20to30_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-20to30_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-20to30_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-15to20_MuEnriched --fileList $dataset_loc/dataset_QCD_Pt-15to20_MuEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-15to20_MuEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

### electron
create-batch --jobName QCD_Pt-300toInf_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-300toInf_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-300toInf_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
#create-batch --jobName QCD_Pt-170to300_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-170to300_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-170to300_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-120to170_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-120to170_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-120to170_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-80to120_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-80to120_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-80to120_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-50to80_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-50to80_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-50to80_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-30to50_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-30to50_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-30to50_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-20to30_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-20to30_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-20to30_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
create-batch --jobName QCD_Pt-15to20_EMEnriched --fileList $dataset_loc/dataset_QCD_Pt-15to20_EMEnriched.txt --cfg $cfg --transferDest $save_loc/QCD_Pt-15to20_EMEnriched --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'


# Systematics
### hdamp up
create-batch --jobName TTLJ_PowhegPythiai_SYS_hdampUp_ttbb --fileList $dataset_loc/dataset_TTLJ_powheg_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampUp_ttbb --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampUp_ttbj --fileList $dataset_loc/dataset_TTLJ_powheg_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampUp_ttbj --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampUp_ttcc --fileList $dataset_loc/dataset_TTLJ_powheg_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampUp_ttcc --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampUp_ttLF --fileList $dataset_loc/dataset_TTLJ_powheg_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampUp_ttLF --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampUp_ttother --fileList $dataset_loc/dataset_TTLJ_powheg_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampUp_ttother --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=5'
create-batch --jobName TTLJ_PowhegPythiaBkg_SYS_hdampUp --fileList $dataset_loc/dataset_TTLJ_powheg_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythiaBkg_SYS_hdampUp --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTLL_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTTo2L2Nu_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTLL_PowhegPythiaBkg_SYS_hdampUp --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTJJ_PowhegPythiaBkg_SYS_hdampUp --fileList $dataset_loc/dataset_TTToHadronic_hdampUP.txt --cfg $cfg --transferDest $save_loc/TTJJ_PowhegPythiaBkg_SYS_hdampUp --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
### hdamp down
create-batch --jobName TTLJ_PowhegPythiai_SYS_hdampDown_ttbb --fileList $dataset_loc/dataset_TTLJ_powheg_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampDown_ttbb --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampDown_ttbj --fileList $dataset_loc/dataset_TTLJ_powheg_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampDown_ttbj --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampDown_ttcc --fileList $dataset_loc/dataset_TTLJ_powheg_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampDown_ttcc --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampDown_ttLF --fileList $dataset_loc/dataset_TTLJ_powheg_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampDown_ttLF --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'
create-batch --jobName TTLJ_PowhegPythia_SYS_hdampDown_ttother --fileList $dataset_loc/dataset_TTLJ_powheg_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_hdampDown_ttother --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=5'
create-batch --jobName TTLJ_PowhegPythiaBkg_SYS_hdampDown --fileList $dataset_loc/dataset_TTLJ_powheg_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythiaBkg_SYS_hdampDown --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTLL_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTTo2L2Nu_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTLL_PowhegPythiaBkg_SYS_hdampDown --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTJJ_PowhegPythiaBkg_SYS_hdampDown --fileList $dataset_loc/dataset_TTToHadronic_hdampDOWN.txt --cfg $cfg --transferDest $save_loc/TTJJ_PowhegPythiaBkg_SYS_hdampDown --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
### TuneCP5 up
create-batch --jobName TTLJ_PowhegPythiai_SYS_TuneCP5Up_ttbb --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Up_ttbb --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Up_ttbj --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Up_ttbj --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Up_ttcc --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Up_ttcc --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Up_ttLF --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Up_ttLF --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Up_ttother --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Up_ttother --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=5'
create-batch --jobName TTLJ_PowhegPythiaBkg_SYS_TuneCP5Up --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythiaBkg_SYS_TuneCP5Up --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTLL_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTTo2L2Nu_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTLL_PowhegPythiaBkg_SYS_TuneCP5Up --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTJJ_PowhegPythiaBkg_SYS_TuneCP5Up --fileList $dataset_loc/dataset_TTToHadronic_TuneCP5up.txt --cfg $cfg --transferDest $save_loc/TTJJ_PowhegPythiaBkg_SYS_TuneCP5Up --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
### TuneCP5 down
create-batch --jobName TTLJ_PowhegPythiai_SYS_TuneCP5Down_ttbb --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Down_ttbb --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Down_ttbj --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Down_ttbj --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Down_ttcc --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Down_ttcc --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Down_ttLF --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Down_ttLF --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'
create-batch --jobName TTLJ_PowhegPythia_SYS_TuneCP5Down_ttother --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythia_SYS_TuneCP5Down_ttother --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=5'
create-batch --jobName TTLJ_PowhegPythiaBkg_SYS_TuneCP5Down --fileList $dataset_loc/dataset_TTLJ_powheg_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLJ_PowhegPythiaBkg_SYS_TuneCP5Down --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTLL_PowhegPythiaBkg --fileList $dataset_loc/dataset_TTTo2L2Nu_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTLL_PowhegPythiaBkg_SYS_TuneCP5Down --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
create-batch --jobName TTJJ_PowhegPythiaBkg_SYS_TuneCP5Down --fileList $dataset_loc/dataset_TTToHadronic_TuneCP5down.txt --cfg $cfg --transferDest $save_loc/TTJJ_PowhegPythiaBkg_SYS_TuneCP5Down --maxFiles 50 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0'
