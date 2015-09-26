#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py 'maxEvents=100' 'runOnRelVal=1' || die 'Failure to run PAT2CAT_cfg.py' $?
