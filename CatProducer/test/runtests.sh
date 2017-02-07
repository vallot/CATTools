#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

INPUTFILE=`python <<EOF
from CATTools.Validation.commonTestInput_cff import *
print commonTestMiniAODs["bkg"][0]
EOF
`
cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6' \
  "inputFiles=$INPUTFILE" \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?

INPUTFILE=`python <<EOF
from CATTools.Validation.commonTestInput_cff import *
print commonTestMiniAODs["sig"][0]
EOF
`
cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6' \
  "inputFiles=$INPUTFILE" \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?

INPUTFILE=`python <<EOF
from CATTools.Validation.commonTestInput_cff import *
print commonTestMiniAODs["data"][0]
EOF
`
cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=0' 'useMiniAOD=1' 'globalTag=80X_dataRun2_2016SeptRepro_v4' \
  "inputFiles=$INPUTFILE" \
  || die 'Failure to run PAT2CAT from MiniAOD real data' $?

