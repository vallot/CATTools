#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_6_2/RelValZMM_13/MINIAODSIM/76X_mcRun2_asymptotic_v12-v1/00000/027AB257-A09C-E511-8F3D-003048FFD736.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_6_2/RelValTTbar_13/MINIAODSIM/76X_mcRun2_asymptotic_v12-v1/00000/0A10678B-8D9C-E511-8AB3-0CC47A4D75F6.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=0' 'useMiniAOD=1' \
  'inputFiles=/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/22890CFA-034A-E511-A23D-02163E011EC1.root' \
  || die 'Failure to run PAT2CAT from MiniAOD real data' $?

