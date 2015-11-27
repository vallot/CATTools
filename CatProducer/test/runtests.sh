#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_15/RelValZMM_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/10FF6E32-3C72-E511-87AD-0025905A60B4.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_15/RelValTTbar_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/0253820F-4772-E511-ADD3-002618943856.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=0' 'useMiniAOD=1' \
  'inputFiles=/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/790/00000/22890CFA-034A-E511-A23D-02163E011EC1.root' \
  || die 'Failure to run PAT2CAT from MiniAOD real data' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=0' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_15/RelValZMM_13/GEN-SIM-RECO/PU25ns_74X_mcRun2_asymptotic_v2-v1/00000/18B13146-3872-E511-A382-00261894394D.root' \
  || die 'Failure to run PAT2CAT from AODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=0' 'runGenTop=1' 'runOnRelVal=1' \
  'inputFiles=/store/relval/CMSSW_7_4_12/RelValTTbar_13/GEN-SIM-RECO/PU25ns_74X_mcRun2_asymptotic_v2_v2-v1/00000/006F3660-4B5E-E511-B8FD-0025905B8596.root' \
  || die 'Failure to run PAT2CAT from AODSIM + GenTop' $?
