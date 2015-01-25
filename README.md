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
git checkout tags/cat53x2
cd ..
scram b -j 8

cmsRun $SRT_CMSSW_BASE_SCRAMRTDEL/src/CATTools/CatProducer/prod/crab/runCatupling.py

```

# Submit CRAB jobs using CATTools
## 1. This is based on crab3
 - no need to change configurations since jobs are stored all in remote servers
 - make sure crab3 is setup first
```bash
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

## 2. Using submitCrab3.py
- This script can be used to pass CRAB configuration parameters from the command line depending on the Dataset lists are located on MC/ or RD/ directories for data type.

### 2-1. Usage 
```bash
submitCrab3.py -n <requestName> -i <inputFile> -s
```
To submit jobs add in '-s', without -s, just the job submission command is displayed
### 2-2. Example
- If we want to submit cattuple ttbar and diboson datasets.
```bash
cd $SRT_CMSSW_BASE_SCRAMRTDEL/src/CATTools/CatProducer/prod/crab
submitCrab3.py -i MC/ttbar_dilepton.txt -n catTooltest -s 
```

### 2-3.Submitting jobs
- because the jobs only saved in local storage (where the job was done) we need to get the output manually. 
- this has the benefit of crab jobs not failing due to transfer errors
- to get the output use "crab out -t <taskdir>" or the script below
```bash
python getcrabOut.py
```
