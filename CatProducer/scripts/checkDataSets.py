#!/usr/bin/env python
import os,sys,getopt,shutil

datasets = []

try:
    opts, args = getopt.getopt(sys.argv[1:],"h:i:",["inputFile"])
except getopt.GetoptError:          
    print 'Usage : ./checkDataSets.py -i <inputFile>'
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./checkDataSets.py -i <inputFile>'
        sys.exit()
    elif opt in ("-i", "--inputFile"):
        inputFile = arg
        if os.path.isfile(inputFile):
            lines = open(inputFile)
            datasets = lines.readlines()
        else:
            datasets.append(inputFile)
            
#print datasets
for dataset in datasets:
    if len(dataset) < 10:
        continue
    if dataset.startswith("#"):
        continue
    #print dataset
    out=os.popen("cli --query='%s' --limit=0"%(dataset)).read()
    
    #print out
    if out != dataset:
        print "DATASET not found", dataset,"similar types are"
        splitDatasetName = dataset.strip().split("/")
        dataPhysics = splitDatasetName[1]
        dataType = splitDatasetName[2]
        dataFormat = splitDatasetName[3]
        dataType = dataType[0:-20] +"*"
        os.system("cli --query='/%s/*/%s' --limit=0"%(dataPhysics,dataFormat))
        print
