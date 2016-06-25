#!/bin/sh

function die { echo $1: status $2 ;  exit $2; }

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runOnRelVal=0' \
  'inputFiles=/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/042C33E6-93CC-E511-B7D6-0CC47A4D7662.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM' $?
#  'inputFiles=/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/06532BBC-05C8-E511-A60A-F46D043B3CE5.root' \

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=1' 'useMiniAOD=1' 'runGenTop=1' 'runOnRelVal=0' \
  'inputFiles=/store/mc/RunIIFall15MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/00000/0024B479-17DC-E511-B25F-0025904CF758.root' \
  || die 'Failure to run PAT2CAT from MiniAODSIM + GenTop' $?
#  'inputFiles=/store/mc/RunIIFall15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/00547C97-2FCC-E511-8D75-002590DB91D2.root' \

cmsRun ${LOCAL_TEST_DIR}/../prod/PAT2CAT_cfg.py \
  'maxEvents=100' 'runOnMC=0' 'useMiniAOD=1' \
  'inputFiles=/store/data/Run2015C_25ns/DoubleMuon/MINIAOD/16Dec2015-v1/20000/081A3AE2-ABB5-E511-9A0D-7845C4FC368C.root' \
  || die 'Failure to run PAT2CAT from MiniAOD real data' $?

#  'inputFiles=/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v3/000/258/158/00000/F4B20DC1-366C-E511-90F2-02163E011985.root' \
