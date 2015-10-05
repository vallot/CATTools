#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/MINIAODSIM/74X_mcRun2_asymptotic_realisticBS_v0_2015Jul24PU-v1/00000/487C927D-6732-E511-9FB3-0025905B85AA.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/MINIAODSIM/MCRUN2_74_V9-v1/00000/2403409D-1225-E511-B64E-0025905A6132.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=0' 'useMiniAOD=1' \
  'inputFiles=/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/22890CFA-034A-E511-A23D-02163E011EC1.root' \
  || die 'Failure to run PAT2CAT from MiniAOD real data' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=0' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_6_patch6/RelValZMM_13/GEN-SIM-RECO/74X_mcRun2_asymptotic_realisticBS_v0_2015Jul24PU-v1/00000/3E2885EA-6032-E511-B000-0025905A60B6.root' \
  || die 'Failure to run PAT2CAT from AODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=0' 'runGenTop=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_6_patch6/RelValTTbar_13/GEN-SIM-RECO/MCRUN2_74_V9-v1/00000/54F6E09C-1225-E511-842B-0025905A612E.root' \
  || die 'Failure to run PAT2CAT from AODSIM + GenTop' $?
