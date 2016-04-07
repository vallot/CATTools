#!/bin/bash
hostname
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 project -n fit CMSSW CMSSW_7_6_3_patch2
cd fit
pwd
eval `scramv1 runtime -sh`
cd ..
./op.py $1 $2
rm -rf fit
