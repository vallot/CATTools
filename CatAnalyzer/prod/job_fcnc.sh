create-batch --jobName SingleMuon_Run2017B  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017B.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleMuon_Run2017B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleMuon_Run2017C  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017C.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleMuon_Run2017C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleMuon_Run2017D  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017D.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleMuon_Run2017D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleMuon_Run2017E  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017E.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleMuon_Run2017E --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleMuon_Run2017F  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017F.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleMuon_Run2017F --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleElectron_Run2017B  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017B.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleElectron_Run2017B --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleElectron_Run2017C  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017C.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleElectron_Run2017C --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleElectron_Run2017D  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017D.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleElectron_Run2017D --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleElectron_Run2017E  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017E.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleElectron_Run2017E --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleElectron_Run2017F  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017F.txt  --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleElectron_Run2017F --maxFiles 30 --args 'UserJSON=true runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName TT_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_powheg_ttbb  --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=1 Powheg=true'

create-batch --jobName TT_powheg_ttbj --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_powheg_ttbj --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=2 Powheg=true'

create-batch --jobName TT_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_powheg_ttcc --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=3 Powheg=true'

create-batch --jobName TT_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_powheg_ttlf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=4 Powheg=true'

create-batch --jobName TT_powheg_ttother --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_powheg_ttother --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=5 Powheg=true'

create-batch --jobName TTLL_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_PSWeight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TTLL_powheg  --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=7 Powheg=true'

create-batch --jobName TTHad_powheg --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_PSWeight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TTHad_powheg --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=7 Powheg=true'

create-batch --jobName DYJets --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_part2.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_part2 --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_4to50_HT70to100 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_4to50_HT70to100.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_4to50_HT70to100 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_4to50_HT100to200 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_4to50_HT100to200.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_4to50_HT100to200 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_4to50_HT100to200_part2 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_4to50_HT100to200_part2.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_4to50_HT100to200_part2 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_4to50_HT200to400 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_4to50_HT200to400.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_4to50_HT200to400 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_4to50_HT400to600 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_4to50_HT400to600.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_4to50_HT400to600 --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName DYJets_4to50_HT600toinf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_4to50_HT600toinf.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/DYJets_4to50_HT600toinf --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName W1JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/W1JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName W2JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/W2JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName W3JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W3JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/W3JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName W4JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/W4JetsToLNu --maxFiles 30 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName ZZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/ZZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName WZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/WZ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName WW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WW.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/WW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleTop_s --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_s.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleTop_s --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=false'

create-batch --jobName SingleTop_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_t.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleTop_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=true'

create-batch --jobName SingleTbar_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_t.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleTbar_t --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=true'

create-batch --jobName SingleTop_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_tW.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleTop_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=true'

create-batch --jobName SingleTbar_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_tW.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/SingleTbar_tW --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=0 TTbarCatMC=0 Powheg=true'

create-batch --jobName TTWJetsToLNu_PSweight --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToLNu_PSweight.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TTWJetsToLNu_PSweight --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0 Powheg=false'

create-batch --jobName TTWJetsToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTWJetsToQQ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TTWJetsToQQ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0 Powheg=false'

create-batch --jobName TTZToLLNuNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToLLNuNu.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TTZToLLNuNu --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0 Powheg=false'

create-batch --jobName TTZToQQ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTZToQQ.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TTZToQQ --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0 Powheg=false'

create-batch --jobName ttHTobb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHTobb.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/ttHTobb --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0 Powheg=true'

create-batch --jobName ttHToNonbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ttHToNonbb.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/ttHToNonbb --maxFiles 20 --args 'UserJSON=false runOnTTbarMC=2 TTbarCatMC=0 Powheg=true'

create-batch --jobName ST_TH_1L3B_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hct.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/ST_TH_1L3B_Hct --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=7'

create-batch --jobName ST_TH_1L3B_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hut.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/ST_TH_1L3B_Hut --maxFiles 10 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=7'

create-batch --jobName TT_TH_1L3B_aTLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_aTLep_Hut.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_TH_1L3B_aTLep_Hut --maxFiles 5 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=7'

create-batch --jobName TT_TH_1L3B_TLep_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TT_TH_1L3B_TLep_Hut.txt --cfg fcncAnalyzer_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/180922/production/TT_TH_1L3B_TLep_Hut --maxFiles 8 --args 'UserJSON=false runOnTTbarMC=1 TTbarCatMC=7'
