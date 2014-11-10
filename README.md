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
cd CATTools
git checkout tags/cat_5_3_22_5
cd ..
scram b -j 8

cmsRun $SRT_CMSSW_BASE_SCRAMRTDEL/src/CATTools/CatProducer/test/runCatupling.py

```

# Submit CRAB jobs using CATTools
## 1. Modify crabMC.cfg and crabRD.cfg files.
- Some variables will be modify to support specific case. 
  - Madantory : **pset**, **storage_element**, **user_remote_dir**
  - Option : total_number_of_events, events_per_job, eMail
- **Do not use _return_data_ for multicrab.**
  - Large result files can be corructed due to storage limit(100GB). 

## 2. Using genMultiCRAB.py
- This script can be used in order to generate correct "multicrab.cfg"
- Dataset lists are located on MC/ or RD/ directories for data type.

### 2-1. Usage 
```bash
./genMultiCRAB.py [Dataset1.txt] [Dataset2.txt]
```

### 2-2. Example
- If we want to generate multicrab.cfg file to make ntuple about ttbar and diboson datasets.
```bash
./genMultiCRAB.py MC/ttbar_dilepton.txt MC/diboson.txt

multicrab -create
multicrab -submit
multicrab -status
multicrab -get
```
