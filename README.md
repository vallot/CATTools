CATTools
========

for cms analysis

Test file : catTuple.root can be found in https://www.dropbox.com/s/40tvwebdv6g0x1m/catTuple.root?dl=0
```bash
scram p -n cat CMSSW CMSSW_5_3_22_patch1
cd cat/src
cmsenv
git-cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd $SRT_CMSSW_BASE_SCRAMRTDEL/src
git clone git@github.com:vallot/CATTools.git
git checkout tags/CMSSW_5_3_22_1
scram b -j 8

cmsRun $SRT_CMSSW_BASE_SCRAMRTDEL/src/CATTools/CatProducer/test/runCatupling.py

```
