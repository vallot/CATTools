create-batch --jobName SingleMuon_Run2016B --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016B.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2016C --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016C.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2016D --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016D.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2016E --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016E.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016E --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2016F --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016F.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016F --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2016G --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016G.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016G --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleMuon_Run2016H --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2016H.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleMuon_Run2016H --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016B --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016B.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016C --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016C.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016D --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016D.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016E --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016E.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016E --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016F --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016F.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016F --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016G --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016G.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016G --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

create-batch --jobName SingleElectron_Run2016H --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2016H.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleElectron_Run2016H --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0'

#Old Tune
#TT
create-batch --jobName TT_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTTo2L2Nu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTTo2L2Nu  --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLJ_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLJ_powheg --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLJ_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_ttbb.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLJ_ttbb --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_aMC --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_aMC.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_aMC --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_evtgen --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_evtgen.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_evtgen --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_herwig_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_herwig_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_herwig_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_herwig_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_herwig_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_herwig_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#fsr/isr
create-batch --jobName TT_powheg_fsrdown_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_fsrdown_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_fsrdown_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_fsrdown_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_fsrdown_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_fsrdown_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_fsrup_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_fsrup_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_fsrup_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_fsrup_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_fsrup_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_fsrup_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_isrdown_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_isrdown_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_isrdown_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_isrdown_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_isrdown_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_isrdown_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_isrup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_isrup.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_isrup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#tune
create-batch --jobName TT_powheg_tunedown_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_tunedown_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_tunedown_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_tunedown_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_tunedown_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_tunedown_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_tuneup_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_tuneup_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_tuneup_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_tuneup_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_tuneup_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_tuneup_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#hdamp
create-batch --jobName TT_powheg_hdampdown_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_hdampdown_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_hdampdown_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_hdampdown_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_hdampdown_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_hdampdown_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_hdampup_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_hdampup_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_hdampup_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_hdampup_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_powheg_hdampup_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_hdampup_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#CP5
#TT
create-batch --jobName TT_powheg_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_CP5 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5  --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTHad_powheg_CP5 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#TT hdamp
create-batch --jobName TT_powheg_CP5_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_CP5_hdampUP.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_CP5_hdampup --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_hdampup_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_hdampUP_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_hdampup_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_hdampup_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_hdampUP_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_hdampup_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_CP5_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_CP5_hdampUP.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTHad_powheg_CP5_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_CP5_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_CP5_hdampDOWN.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_CP5_hdampdown --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_hdampdown_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_hdampDOWN_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_hdampdown_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_hdampdown_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_hdampDOWN_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_hdampdown_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_CP5_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_CP5_hdampDOWN.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTHad_powheg_CP5_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

#TT tune
create-batch --jobName TT_powheg_CP5_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_CP5_TuneCP5up.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_CP5_TuneCP5up --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_TuneCP5up_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_TuneCP5up_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_TuneCP5up_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_TuneCP5up_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_TuneCP5up_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_TuneCP5up_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_CP5_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_CP5_TuneCP5up.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTHad_powheg_CP5_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TT_powheg_CP5_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_CP5_TuneCP5down.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_powheg_CP5_TuneCP5down --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_TuneCP5down_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_TuneCP5down_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_TuneCP5down_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLL_powheg_CP5_TuneCP5down_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_CP5_TuneCP5down_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLL_powheg_CP5_TuneCP5down_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTHad_powheg_CP5_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_CP5_TuneCP5down.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTHad_powheg_CP5_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1'

create-batch --jobName TTLJ_ttbb_powheg_openloops_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_ttbb_powheg_openloops_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTLJ_ttbb_powheg_openloops_CP5 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=2'

#Other SM
create-batch --jobName DYJets_10to50_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_10to50_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/DYJets_10to50_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName DYJets_10to50_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_10to50_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/DYJets_10to50_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName DYJets --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/DYJets --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W1JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W1JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W2JetsToLNu_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W2JetsToLNu_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W2JetsToLNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W2JetsToLNu_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W3JetsToLNu_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W3JetsToLNu_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W3JetsToLNu_part1 --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W3JetsToLNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W3JetsToLNu_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W3JetsToLNu_part2 --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W4JetsToLNu_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W4JetsToLNu_part1 --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W4JetsToLNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W4JetsToLNu_part2 --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName W4JetsToLNu_part3 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu_part3.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/W4JetsToLNu_part3 --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WJetsToLNu_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WJetsToLNu_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WJetsToLNu_part1 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WJetsToLNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WJetsToLNu_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WJetsToLNu_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_s --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_s.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTop_s --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_t.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTop_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTbar_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_t.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTbar_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_tW.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTop_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTbar_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_tW.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTbar_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_s_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_s_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTop_s_CP5 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_t_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_t_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTop_t_CP5 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTbar_t_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_t_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTbar_t_CP5 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTop_tW_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_tW_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTop_tW_CP5 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName SingleTbar_tW_CP5 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_tW_CP5.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/SingleTbar_tW_CP5 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WW_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WW_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WW_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WW_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WW_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WW_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WZ_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZ_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WZ_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WZ_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZ_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WZ_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName ZZ_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZ_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ZZ_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName ZZ_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZ_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ZZ_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WWW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WWW.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WWW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WWZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WWZ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WWZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName WZZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZZ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/WZZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName ZZZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZZ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ZZZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName TTWJetsToLNu_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToLNu_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTWJetsToLNu_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTWJetsToLNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToLNu_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTWJetsToLNu_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTWJetsToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToQQ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTWJetsToQQ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTZToLLNuNu_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToLLNuNu_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTZToLLNuNu_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTZToLLNuNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToLLNuNu_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTZToLLNuNu_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTZToLLNuNu_part3 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToLLNuNu_part3.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTZToLLNuNu_part3 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTZToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToQQ.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTZToQQ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName ttHTobb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHTobb.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ttHTobb --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName ttHToNonbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHToNonbb.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ttHToNonbb --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName ST_TH_1L3B_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hct.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ST_TH_1L3B_Hct --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName ST_TH_1L3B_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hut.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/ST_TH_1L3B_Hut --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TT_TH_1L3B_aTLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hut.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_TH_1L3B_aTLep_Hut --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TT_TH_1L3B_TLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hut.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_TH_1L3B_TLep_Hut --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TT_TH_1L3B_aTLep_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hct.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_TH_1L3B_aTLep_Hct --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TT_TH_1L3B_TLep_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hct.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TT_TH_1L3B_TLep_Hct --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName TTTT --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTT.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/TTTT --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=2'

create-batch --jobName QCD_EM20to30 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-20to30_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM20to30 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM30to50_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_EMEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM30to50_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM30to50_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_EMEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM30to50_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM50to80_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_EMEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM50to80_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM50to80_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_EMEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM50to80_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM80to120_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_EMEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM80to120_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM80to120_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_EMEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM80to120_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM120to170_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_EMEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM120to170_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM120to170_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_EMEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM120to170_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM170to300 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM170to300 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_EM300toInf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300toInf_EMEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_EM300toInf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu15to20 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-15to20_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu15to20 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu20to30 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-20to30_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu20to30 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu30to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu30to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu50to80 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu50to80 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu80to120_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu80to120_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu80to120_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu80to120_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu120to170 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_MuEnriched.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu120to170 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu170to300_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu170to300_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu170to300_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu170to300_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu300to470_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300to470_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu300to470_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu300to470_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300to470_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu300to470_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu300to470_part3 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300to470_MuEnriched_part3.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu300to470_part3 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu470to600_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-470to600_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu470to600_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu470to600_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-470to600_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu470to600_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu470to600_part3 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-470to600_MuEnriched_part3.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu470to600_part3 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu600to800_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-600to800_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu600to800_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu600to800_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-600to800_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu600to800_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu800to1000_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-800to1000_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu800to1000_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu800to1000_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-800to1000_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu800to1000_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu800to1000_part3 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-800to1000_MuEnriched_part3.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu800to1000_part3 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu1000toInf_part1 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-1000toInf_MuEnriched_part1.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu1000toInf_part1 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'

create-batch --jobName QCD_Mu1000toInf_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-1000toInf_MuEnriched_part2.txt --cfg testAnalyzer_cfg.py --transferDest /store/user/jipark/cat_test_V8_1/QCD_Mu1000toInf_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0'
