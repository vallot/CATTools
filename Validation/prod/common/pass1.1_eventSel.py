#!/usr/bin/env python
import sys, os

## Load JSON file and categorize datasets
import json
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
js = json.loads(open("%s/dataset.json" % dataDir).read())

outDir = "pass1"
if os.path.exists(outDir):
    print "Directory already exists. remove or rename it."
    sys.exit()
os.mkdir(outDir)
os.system("ln -s ../analyze_mc_cfg.py %s/analyze_mc_cfg.py" % outDir)
os.system("ln -s ../analyze_data_cfg.py %s/analyze_data_cfg.py" % outDir)

fout = open("%s/submit.sh" % outDir, 'w')
cmds = {}
for sample in js:
    jobName = sample['name']
    type = 'mc'
    if sample['type'] == 'Data': type = 'data'

    submitCmd  = "create-batch --cfg analyze_%s_cfg.py --maxFiles 1 " % type
    submitCmd += " --jobName {0} --fileList {1}/dataset_{0}.txt ".format(jobName, dataDir)

    print>>fout, submitCmd
os.system("chmod +x %s/submit.sh" % outDir)
