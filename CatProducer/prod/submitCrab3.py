#!/usr/bin/env python
import os,sys

if len(sys.argv) < 2 :
    print "Usage : ./submitCrab3.py ttbar_mc.txt ttbar_rd.txt ..."
    exit(-1)

datasets =[]
for file in sys.argv[1:] :
    lines = open(file)
    datasets += lines.readlines()

### MC or Data?
isMC = True
type_label = "MC"
requestName = "cat72x2"
publishDataName = "top"
runOnMC = "cmsRun runCatupling.py runOnMC=True"

dataset = datasets[0]
datatype = dataset.strip().split("/")[-1]
if datatype == "AOD" or datatype == "MINIAOD" :
    isMC = False
    type_label = "RD"
    runOnMC = "cmsRun runCatupling.py runOnMC=False"

crabcommand ='crab submit -c crabConfigMC.py'

if isMC :
    print "I guess these datasets are Monte Carlo samples.crab3 job is MC.\n"
else :
    print "I guess these datasets are Real Data samples.crab3 job is RD.\n"
    crabcommand ='crab submit -c crabConfigRD.py'
    
for dataset in datasets :
    dataset = dataset.strip()
    if isMC :
        label = dataset.split("/")[1]
    else :
        label = dataset.split("/")[1]+"_"+dataset.split("/")[2]

    sendjob = crabcommand + " Data.publishDataName='%s' General.requestName='%s_%s' Data.inputDataset='%s'"%(publishDataName,requestName,label,dataset)
    print sendjob
    os.system(sendjob)
