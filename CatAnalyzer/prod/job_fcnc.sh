create-batch --jobName SingleMuon_Run2017B  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017B.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleMuon_Run2017B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleMuon_Run2017C  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017C.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleMuon_Run2017C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleMuon_Run2017D  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017D.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleMuon_Run2017D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleMuon_Run2017E  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017E.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleMuon_Run2017E --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleMuon_Run2017F  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017F.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleMuon_Run2017F --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleElectron_Run2017B  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017B.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleElectron_Run2017B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleElectron_Run2017C  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017C.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleElectron_Run2017C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleElectron_Run2017D  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017D.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleElectron_Run2017D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleElectron_Run2017E  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017E.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleElectron_Run2017E --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleElectron_Run2017F  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017F.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleElectron_Run2017F --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0'

#TT
create-batch --jobName TT_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttbb --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TT_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttcc --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TT_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttlf --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTLL_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttbb  --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTLL_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttcc  --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTLL_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttlf  --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTHad_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttbb --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTHad_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttcc --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTHad_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttlf --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

#TT hdamp
create-batch --jobName TT_powheg_ttbb_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttbb_hdampup --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TT_powheg_ttcc_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttcc_hdampup --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TT_powheg_ttlf_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttlf_hdampup --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTLL_powheg_ttbb_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttbb_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1 PUMap=TTTo2L2Nu_hdampUP'

create-batch --jobName TTLL_powheg_ttcc_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttcc_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2 PUMap=TTTo2L2Nu_hdampUP'

create-batch --jobName TTLL_powheg_ttlf_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttlf_hdampup --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3 PUMap=TTTo2L2Nu_hdampUP'

create-batch --jobName TTHad_powheg_ttbb_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttbb_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTHad_powheg_ttcc_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttcc_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTHad_powheg_ttlf_hdampup --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampUP.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttlf_hdampup --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TT_powheg_ttbb_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttbb_hdampdown --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1 PUMap=TTLJ_powheg_hdampDOWN'

create-batch --jobName TT_powheg_ttcc_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttcc_hdampdown --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2 PUMap=TTLJ_powheg_hdampDOWN'

create-batch --jobName TT_powheg_ttlf_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttlf_hdampdown --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3 PUMap=TTLJ_powheg_hdampDOWN'

create-batch --jobName TTLL_powheg_ttbb_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttbb_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTLL_powheg_ttcc_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttcc_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTLL_powheg_ttlf_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttlf_hdampdown --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTHad_powheg_ttbb_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttbb_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTHad_powheg_ttcc_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttcc_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTHad_powheg_ttlf_hdampdown --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_hdampDOWN.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttlf_hdampdown --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

#TT tune
create-batch --jobName TT_powheg_ttbb_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttbb_TuneCP5up --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TT_powheg_ttcc_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttcc_TuneCP5up --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TT_powheg_ttlf_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttlf_TuneCP5up --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTLL_powheg_ttbb_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttbb_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTLL_powheg_ttcc_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttcc_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTLL_powheg_ttlf_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttlf_TuneCP5up --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTHad_powheg_ttbb_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttbb_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTHad_powheg_ttcc_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttcc_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTHad_powheg_ttlf_TuneCP5up --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5up.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttlf_TuneCP5up --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TT_powheg_ttbb_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttbb_TuneCP5down --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1 PUMap=TTLJ_powheg_TuneCP5down'

create-batch --jobName TT_powheg_ttcc_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttcc_TuneCP5down --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2 PUMap=TTLJ_powheg_TuneCP5down'

create-batch --jobName TT_powheg_ttlf_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_powheg_ttlf_TuneCP5down --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3 PUMap=TTLJ_powheg_TuneCP5down'

create-batch --jobName TTLL_powheg_ttbb_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttbb_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTLL_powheg_ttcc_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttcc_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTLL_powheg_ttlf_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTLL_powheg_ttlf_TuneCP5down --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

create-batch --jobName TTHad_powheg_ttbb_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttbb_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1'

create-batch --jobName TTHad_powheg_ttcc_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttcc_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2'

create-batch --jobName TTHad_powheg_ttlf_TuneCP5down --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_TuneCP5down.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTHad_powheg_ttlf_TuneCP5down --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3'

#Other SM
create-batch --jobName DYJets_10to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_10to50.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/DYJets_10to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=DYJets_10to50'

create-batch --jobName DYJets --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/DYJets --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName W1JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/W1JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName W2JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/W2JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName W3JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W3JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/W3JetsToLNu --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=W3JetsToLNu'

create-batch --jobName W4JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/W4JetsToLNu --maxFiles 15 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName WJetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WJetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/WJetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName WJetsToLNu_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WJetsToLNu_part2.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/WJetsToLNu_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleTop_s --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_s.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleTop_s --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleTop_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_t.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleTop_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleTbar_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_t.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleTbar_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleTop_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_tW.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleTop_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName SingleTbar_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_tW.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/SingleTbar_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName WW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WW.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/WW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=WW'

create-batch --jobName WZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/WZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=WZ'

create-batch --jobName ZZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/ZZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=ZZ'

create-batch --jobName TTWJetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTWJetsToLNu --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName TTWJetsToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToQQ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTWJetsToQQ --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName TTZToLLNuNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToLLNuNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTZToLLNuNu --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4 PUMap=TTZToLLNuNu'

create-batch --jobName TTZToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToQQ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TTZToQQ --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName ttHTobb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHTobb.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/ttHTobb --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName ttHToNonbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHToNonbb.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/ttHToNonbb --maxFiles 2 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName ST_TH_1L3B_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hct.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/ST_TH_1L3B_Hct --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName ST_TH_1L3B_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hut.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/ST_TH_1L3B_Hut --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName TT_TH_1L3B_aTLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hut.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_TH_1L3B_aTLep_Hut --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName TT_TH_1L3B_TLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hut.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_TH_1L3B_TLep_Hut --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName TT_TH_1L3B_aTLep_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hct.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_TH_1L3B_aTLep_Hct --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName TT_TH_1L3B_TLep_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hct.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/TT_TH_1L3B_TLep_Hct --maxFiles 4 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4'

create-batch --jobName QCD_EM15to20 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-15to20_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM15to20 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-15to20_EMEnriched'

create-batch --jobName QCD_EM20to30 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-20to30_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM20to30 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-20to30_EMEnriched'

create-batch --jobName QCD_EM30to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM30to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-30to50_EMEnriched'

create-batch --jobName QCD_EM50to80 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM50to80 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-50to80_EMEnriched'

create-batch --jobName QCD_EM80to120 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM80to120 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_EM120to170 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM120to170 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-300toInf_EMEnriched'

#create-batch --jobName QCD_EM170to300 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM170to300 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_EM300toInf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300toInf_EMEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_EM300toInf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_Mu15to20 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-15to20_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu15to20 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_Mu20to30 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-20to30_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu20to30 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-20to30_MuEnriched'

create-batch --jobName QCD_Mu30to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-30to50_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu30to50 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-30to50_MuEnriched'

create-batch --jobName QCD_Mu50to80 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-50to80_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu50to80 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-50to80_MuEnriched'

create-batch --jobName QCD_Mu80to120 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-80to120_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu80to120 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-80to120_MuEnriched'

create-batch --jobName QCD_Mu120to170 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-120to170_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu120to170 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-120to170_MuEnriched'

create-batch --jobName QCD_Mu170to300 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-170to300_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu170to300 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-170to300_MuEnriched'

create-batch --jobName QCD_Mu300to470 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-300to470_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu300to470 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_Mu470to600 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-470to600_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu470to600 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 PUMap=QCD_Pt-470to600_MuEnriched'

create-batch --jobName QCD_Mu600to800 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-600to800_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu600to800 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_Mu800to1000 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-800to1000_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu800to1000 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'

create-batch --jobName QCD_Mu1000toInf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_QCD_Pt-1000toInf_MuEnriched.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/jipark/ntuple_jw/V9_7/201201/production/QCD_Mu1000toInf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0'
