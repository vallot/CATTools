CATTools
========

for cms analysis

Test file : catTuple.root can be found in https://www.dropbox.com/s/40tvwebdv6g0x1m/catTuple.root?dl=0
```bash
scram p -n cat CMSSW CMSSW_7_2_1_patch2
cd cat/src
cmsenv
git clone git@github.com:vallot/CATTools.git
cd CATTools
git checkout CMSSW_7_2_X
cd ..
scram b -j 8

cmsRun $SRT_CMSSW_BASE_SCRAMRTDEL/src/CATTools/CatProducer/test/runCatupling.py

```
