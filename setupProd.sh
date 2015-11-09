#!/bin/bash

scram p -n cattools CMSSW CMSSW_7_4_15
cd cattools/src
cmsenv
git cms-merge-topic vallot:cat75x
git cms-addpkg RecoMET/METPUSubtraction
git clone https://github.com/rfriese/RecoMET-METPUSubtraction RecoMET/METPUSubtraction/data -b 74X-13TeV-50nsData-September2015
rm -rf RecoMET/METPUSubtraction/data/.git

git clone git@github.com:vallot/CATTools.git
cd CATTools
git submodule init
git submodule update
git checkout -b prod v7-4-5

cd $CMSSW_BASE/src
scram setup lhapdf
scram b -j 20

