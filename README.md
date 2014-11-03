CATTools
========

for cms analysis
Test file : /pnfs/user/jlee/catTuple.root
```bash
scram p -n cat CMSSW CMSSW_5_3_22_patch1
cmsenv
git-cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -
git clone https://github.com/$USER/CATTools.git
git checkout CMSSW_5_3_22
scram b -j 8
```
