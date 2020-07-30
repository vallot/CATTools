create-batch --jobName SingleMuon_Run2018A --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2018A.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleMuon_Run2018A --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2018B --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2018B.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleMuon_Run2018B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2018C --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2018C.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleMuon_Run2018C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2018D --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2018D.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleMuon_Run2018D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2018A --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_EGamma_Run2018A.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleElectron_Run2018A --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2018B --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_EGamma_Run2018B.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleElectron_Run2018B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2018C --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_EGamma_Run2018C.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleElectron_Run2018C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2018D --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_EGamma_Run2018D.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleElectron_Run2018D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

#TT
create-batch --jobName TT_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_powheg --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTLL_powheg  --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTHad_powheg --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#TT hdamp
create-batch --jobName TT_powheg_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampUP.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_powheg_ttbb_hdampup --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampUP.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTLL_powheg_ttbb_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampUP.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTHad_powheg_ttbb_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampDOWN.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_powheg_ttbb_hdampdown --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampDOWN.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTLL_powheg_ttbb_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampDOWN.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTHad_powheg_ttbb_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#TT tune
create-batch --jobName TT_powheg_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5up.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_powheg_ttbb_TuneCP5up --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5up.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTLL_powheg_ttbb_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5up.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTHad_powheg_ttbb_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5down.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_powheg_ttbb_TuneCP5down --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5down.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTLL_powheg_ttbb_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5down.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTHad_powheg_ttbb_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#Other SM
create-batch --jobName DYJets_10to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_10to50.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/DYJets_10to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName DYJets --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/DYJets --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W1JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/W1JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W2JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/W2JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W3JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W3JetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/W3JetsToLNu --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W4JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/W4JetsToLNu --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WJetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WJetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/WJetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_s --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_s.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleTop_s --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_t.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleTop_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTbar_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_t.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleTbar_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_tW.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleTop_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTbar_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_tW.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/SingleTbar_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WW.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/WW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/WZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName ZZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/ZZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName TTWJetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTWJetsToLNu --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTWJetsToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToQQ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTWJetsToQQ --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTZToLLNuNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToLLNuNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTZToLLNuNu --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTZToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToQQ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TTZToQQ --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName ttHTobb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHTobb.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/ttHTobb --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName ttHToNonbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHToNonbb.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/ttHToNonbb --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName ST_TH_1L3B_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hct.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/ST_TH_1L3B_Hct --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName ST_TH_1L3B_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hut.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/ST_TH_1L3B_Hut --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_TH_1L3B_aTLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hut.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_TH_1L3B_aTLep_Hut --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_TH_1L3B_TLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hut.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_TH_1L3B_TLep_Hut --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_TH_1L3B_aTLep_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hct.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_TH_1L3B_aTLep_Hct --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_TH_1L3B_TLep_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hct.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/TT_TH_1L3B_TLep_Hct --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName QCD_EM15to20 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-15to20_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM15to20 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM20to30 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-20to30_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM20to30 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM30to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM30to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM50to80 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM50to80 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM80to120 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM80to120 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM120to170 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM120to170 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM170to300 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM170to300 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM300toInf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300toInf_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_EM300toInf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu15to20 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-15to20_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu15to20 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu20to30 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-20to30_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu20to30 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu30to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu30to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu50to80 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu50to80 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu80to120 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu80to120 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu120to170 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu120to170 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu170to300 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu170to300 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu300to470 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300to470_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu300to470 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu470to600 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-470to600_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu470to600 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu600to800 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-600to800_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu600to800 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu800to1000 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-800to1000_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu800to1000 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu1000toInf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-1000toInf_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V10_3/QCD_Mu1000toInf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'
