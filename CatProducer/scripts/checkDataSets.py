#!/usr/bin/env python
import os,sys,getopt,shutil,json

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json"%os.environ['CMSSW_BASE']))

for d in datasets:
    dataset = d['DataSetName']
    #print "cli --query='%s' --limit=0"%(dataset)
    out=os.popen("cli --query='%s' --limit=0"%(dataset)).read()
    if out.rstrip() != dataset:
        print "DATASET not found", dataset,"similar types are"
        splitDatasetName = dataset.strip().split("/")
        dataPhysics = splitDatasetName[1]
        dataType = splitDatasetName[2]
        dataFormat = splitDatasetName[3]
        dataType = dataType[0:-20] +"*"
        os.system("cli --query='/%s/*/%s' --limit=0"%(dataPhysics,dataFormat))
        print
