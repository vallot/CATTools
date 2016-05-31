CATTools
========

for cms analysis

# Links
 - CATTuple production
  - Campaigns, recipe : https://github.com/vallot/CATTools/wiki
  - Summary table : https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit#gid=513081606

```bash
# See https://www.cms-kr.org/twiki/bin/view/Computing/CATTools/DataSets/WebHome for installation
cd $CMSSW_BASE/src/CATTools/CatProducer/prod

cmsRun PAT2CAT_cfg.py

# or

cmsRun PAT2CAT_cfg.py useMiniAOD=True inputFiles=/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root

# or

cmsRun PAT2CAT_cfg.py useMiniAOD=False inputFiles=/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root globalTag='PHYS14_25_V2::All'
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
