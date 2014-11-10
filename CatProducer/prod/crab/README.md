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
```
