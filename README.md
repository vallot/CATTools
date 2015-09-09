CATTools
========

for cms analysis

Test file : catTuple.root can be found in /afs/cern.ch/user/j/jlee/public/catTuple.root
```bash
scram p -n cat74 CMSSW CMSSW_7_4_7_patch2
cd cat74/src
cmsenv
git cms-addpkg TopQuarkAnalysis/TopEventProducers
git cms-addpkg CommonTools/PileupAlgos
git cms-addpkg RecoMET/METPUSubtraction
git cms-merge-topic jhgoh:PseudoTop
git cms-merge-topic nhanvtran:puppi-etadep-746p2-v8
git cms-merge-topic cms-met:METCorUnc74X
git clone https://github.com/cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_74X
git clone https://github.com/rfriese/RecoMET-METPUSubtraction RecoMET/METPUSubtraction/data -b 74X-13TeV-Summer15-July2015
rm -rf RecoMET/METPUSubtraction/data/.git
git cms-merge-topic ikrav:egm_id_747_v2
git cms-merge-topic YoungKwonJo:NewGenHFHadronMatcher_cmssw747p2

git clone https://github.com/vallot/CATTools.git -n
cd CATTools
git submodule init
git submodule update
git checkout -b v740 v7-4-0
cd ..

scram setup lhapdf
scram b -j 8

cd CATTools/CatProducer/prod

cmsRun PAT2CAT_cfg.py 

or 

cmsRun PAT2CAT_cfg.py useMiniAOD=True inputFiles=/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root

or 

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
