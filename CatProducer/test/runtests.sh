#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runOnRelVal=0' \
  'inputFiles=/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/004A0552-3929-E611-BD44-0025905A48F0.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'runOnRelVal=1' \
  'inputFiles=/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=0' 'useMiniAOD=1' \
  'inputFiles=/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/275/375/00000/FEDDBC5C-4339-E611-8740-02163E0136C4.root' \
  || die 'Failure to run PAT2CAT from MiniAOD real data' $?

