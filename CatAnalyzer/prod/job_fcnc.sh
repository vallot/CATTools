./create-batch --jobName SingleMuon_Run2017B  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017B.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleMuon_Run2017B --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleMuon_Run2017C  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017C.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleMuon_Run2017C --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleMuon_Run2017D  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017D.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleMuon_Run2017D --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleMuon_Run2017E  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017E.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleMuon_Run2017E --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleMuon_Run2017F  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleMuon_Run2017F.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleMuon_Run2017F --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleElectron_Run2017B  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017B.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleElectron_Run2017B --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleElectron_Run2017C  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017C.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleElectron_Run2017C --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleElectron_Run2017D  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017D.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleElectron_Run2017D --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleElectron_Run2017E  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017E.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleElectron_Run2017E --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleElectron_Run2017F  --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleElectron_Run2017F.txt  --cfg fcncAnalyzer_RD_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleElectron_Run2017F --maxFiles 30 --args 'UserJSON=true' 'runOnTTbarMC=0' 'TTbarCatMC=0'

/create-batch --jobName TT_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_ttbb_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TT_powheg_ttbb  --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=1'

./create-batch --jobName TT_powheg_ttbj --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_ttbj_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TT_powheg_ttbj --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=2'

./create-batch --jobName TT_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_ttcc_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TT_powheg_ttcc --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=3'

./create-batch --jobName TT_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcnvAnalyzer_ttlf_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TT_powheg_ttlf --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=4'

./create-batch --jobName TT_powheg_ttother --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTLJ_powheg_PSweight.txt --cfg fcncAnalyzer_ttother_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/TT_powheg_ttother --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=5'

./create-batch --jobName TTLL_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_PSWeight.txt --cfg fcncAnalyzer_ttbb_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTLL_powheg_ttbb  --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=1'

./create-batch --jobName TTLL_powheg_ttbj --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_PSWeight.txt --cfg fcncAnalyzer_ttbj_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTLL_powheg_ttbj --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=2'

./create-batch --jobName TTLL_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_PSWeight.txt --cfg fcncAnalyzer_ttcc_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTLL_powheg_ttcc --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=3'

./create-batch --jobName TTLL_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_PSWeight.txt --cfg fcnvAnalyzer_ttlf_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTLL_powheg_ttlf --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=4'

./create-batch --jobName TTLL_powheg_ttother --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTTo2L2Nu_PSWeight.txt --cfg fcncAnalyzer_ttother_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTLL_powheg_ttother --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=5'

./create-batch --jobName TTHad_powheg_ttbb --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_PSWeight.txt --cfg fcncAnalyzer_ttbb_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTHad_powheg_ttbb  --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=1'

./create-batch --jobName TTHad_powheg_ttbj --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_PSWeight.txt --cfg fcncAnalyzer_ttbj_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTHad_powheg_ttbj --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=2'

./create-batch --jobName TTHad_powheg_ttcc --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_PSWeight.txt --cfg fcncAnalyzer_ttcc_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTHad_powheg_ttcc --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=3'

./create-batch --jobName TTHad_powheg_ttlf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_PSWeight.txt --cfg fcnvAnalyzer_ttlf_cfg.py --transferDest /store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTHad_powheg_ttlf --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=4'

./create-batch --jobName TTHad_powheg_ttother --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_TTToHadronic_PSWeight.txt --cfg fcncAnalyzer_ttother_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/TTHad_powheg_ttother --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=5'

./create-batch --jobName DYJets --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/DYJets --maxFiles 40 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName DYJets_10to50 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_DYJets_10to50.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/DYJets_10to50 --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W1JetsToLNu_50-150 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu_50-150.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W1JetsToLNu_50-150 --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W1JetsToLNu_150-250 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu_150-250.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W1JetsToLNu_150-250 --maxFiles 30 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W1JetsToLNu_250-400 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu_250-400.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W1JetsToLNu_250-400 --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W1JetsToLNu_400-inf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W1JetsToLNu_400-inf.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W1JetsToLNu_400-inf --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W2JetsToLNu_250-400 --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu_250-400.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W1JetsToLNu_250-400 --maxFiles 30 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W2JetsToLNu_400-inf --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W2JetsToLNu_400-inf.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W2JetsToLNu_400-inf --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W3JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W3JetsToLNu.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W3JetsToLNu --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName W4JetsToLNu --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_W4JetsToLNu.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/W4JetsToLNu --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName ZZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ZZ.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/ZZ --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName WZ --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WZ.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/WZ --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName WW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_WW.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/WW --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleTop_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_t.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleTop_t --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleTbar_t --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_t.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleTbar_t --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleTop_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTop_tW.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleTop_tW --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName SingleTbar_tW --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_SingleTbar_tW.txt --cfg ttbbLepJetsAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/2017/V9_2/production/SingleTbar_tW --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=0' 'TTbarCatMC=0'

./create-batch --jobName ST_TH_1L3B_Hct --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hct.txt --cfg fcncAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/ --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=0'

./create-batch --jobName ST_TH_1L3B_Hut --fileList $CMSSW_BASE/src/CATTools/CatAnalyzer/data/dataset/dataset_ST_TH_1L3B_Hut.txt --cfg fcncAnalyzer_MC_cfg.py --transferDest /xrootd/store/user/minerva1993/ntuple_jw/ --maxFiles 20 --args 'UserJSON=false' 'runOnTTbarMC=1' 'TTbarCatMC=0'
