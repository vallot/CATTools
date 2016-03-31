#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

analysis = 'h2muAnalyzer'
#analysis = 'TtbarDiLeptonAnalyzer'

pythonCfg = 'run_'+analysis+'_cfg.py'
#analysis=analysis+'Silver'

import os,json
datadir = os.environ["CMSSW_BASE"]+'/src/CATTools/CatAnalyzer/data/dataset/'
dataset_json = datadir + 'dataset.json'

with open(dataset_json) as data_file:    
    data = json.load(data_file)
    for i in data:
        #print data[0]
        datasetName = i['name']
        if "QCD" in datasetName:
            continue
        if "ttH" in datasetName:
            continue       
        if "ttZ" in datasetName:
            continue
        if "ttW" in datasetName:
            continue

        fileList = datadir + 'dataset_' + datasetName + '.txt'
        jobName = analysis+'_'+datasetName
        createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 10"%(pythonCfg, jobName, fileList)
        print createbatch
        os.system(createbatch)
        
