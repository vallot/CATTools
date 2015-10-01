#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runOnRelVal=1' \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'runOnRelVal=1' \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=0' 'runOnRelVal=1' \
  || die 'Failure to run PAT2CAT from AODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=0' 'runGenTop=1' 'runOnRelVal=1' \
  || die 'Failure to run PAT2CAT from AODSIM + GenTop' $?
