CATTools
========

for cms analysis

```bash
scram p -n cat CMSSW CMSSW_7_0_6
cd cat/src
cmsenv
git clone git@github.com:vallot/CATTools.git
cd CATTools
git checkout 70X_v1
cd ..
scram b -j 8

cmsRun $SRT_CMSSW_BASE_SCRAMRTDEL/src/CATTools/CatProducer/prod/MiniAOD2CAT_cfg.py

```
