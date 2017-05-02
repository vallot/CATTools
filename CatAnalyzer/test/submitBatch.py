#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

analysis = 'TtbarDiLeptonAnalyzer'
import os
username = os.environ['USER']
pythonCfg = 'run_'+analysis+'_cfg.py'
#analysis=analysis+'Silver'

import os,json
datadir = os.environ["CMSSW_BASE"]+'/src/CATTools/CatAnalyzer/data/dataset/'
dataset_json = datadir + 'dataset.json'
version = os.environ["CMSSW_VERSION"]
mclist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
datalist = ['DoubleEG','DoubleMuon','MuonEG']


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
        #if not "SingleMuon" in datasetName:
        #    continue
        #if not any([d in datasetName for d in datalist]):
        #            continue
        if not any([d in datasetName for d in mclist]):
                    continue
        fileList = datadir + 'dataset_' + datasetName + '.txt'
        jobName = analysis+'_'+datasetName
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 10"%(pythonCfg, jobName, fileList) 
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 50 --transferDest \"root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/pseudotop/ntuples/%s/%s\""%(pythonCfg, jobName, fileList, version, datasetName)
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 50 --transferDest \"root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/pseudotop/ntuples/%s_bump2/%s\""%(pythonCfg, jobName, fileList, version, datasetName)
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 20 --transferDest \"root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/king11kr/ntuples_TtbarDstar_v806/%s_v806/%s\""%(pythonCfg, jobName, fileList, version, datasetName) 
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 20 "%(pythonCfg, jobName, fileList, version, datasetName) 
        createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 20 --transferDest \"root://cms-xrdr.sdfarm.kr:1094//xrd/store/user/%s/ntuples_TtbarDstar_v806/%s/%s\""%(pythonCfg, jobName, fileList, username, version, datasetName)
        print createbatch
        #os.system(createbatch)
        
