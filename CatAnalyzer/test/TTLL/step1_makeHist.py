#!/usr/bin/env python

import sys, os
outDir = "step1_makeHist"
if os.path.exists(outDir):
    print "Directory already exists. remove or rename it."
    os.exit()
os.mkdir(outDir)

## Load JSON file and categorize datasets
import json
dataDir = "%s/src/CATTools/CatAnalyzer/data" % os.environ["CMSSW_BASE"]
js = json.loads(open("%s/dataset.json" % dataDir).read())

sigList = []
bkgList = []
dataList = []

for j in js:
    datasetName = j['DataSetName']
    name = j['name']
    title = j['title']
    if datasetName.endswith("/MINIAOD"): dataList.append(j)
    elif title.startswith('t#bar{t}'): sigList.append(j)
    else: bkgList.append(j)

## Undertainty variations
systAll = {
    'lep_pt/up':'muon.scaleDirection=1 electron.scaleDirection=1',
    'lep_pt/dn':'muon.scaleDirection=-1 electron.scaleDirection=-1',
    'jet_cor/up':'jet.scaleDirection=1',
    'jet_cor/dn':'jet.scaleDirection=-1',
}

systMC = {
    'jet_res/up':'jet.resolDirection=1',
    'jet_res/dn':'jet.resolDirection=-1',
    'pileup/up':'vertex.pileupWeight="pileupWeight:up"',
    'pileup/dn':'vertex.pileupWeight="pileupWeight:dn"',
}

systSig_aMC = {}
for i in range(1, 11):
    systSig_aMC['gen_scale/%d' % i] = 'src="genWeight:pdfWeights" genWeight.index=%d' % i
for i in range(11, 112):
    systSig_aMC['gen_PDF/%d' % i] = ('src="genWeight:pdfWeights" genWeight.index=%d' % i)

## Write script to run create-batch
out = open("%s/submit.sh" % outDir, "w")
for d in sigList:
    name = d['name']
    title = d['title']
    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 100"
    submitCmd += " --fileList %s/dataset_%s.txt" % (dataDir, name)

    print>>out, submitCmd
    ## Loop over all systematics
    for systName in systAll:
        arg = systAll[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
    for systName in systMC:
        arg = systMC[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
    if 'aMC' in title:
        for systName in systSig_aMC:
            arg = systSig_aMC[systName]
            print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))

for d in bkgList:
    name = d['name']
    submitCmd  = "create-batch --cfg analyze_bkg_cfg.py --maxFiles 100"
    submitCmd += " --jobName %s --fileList %s/dataset_%s.txt" % (name, dataDir, name)

    print>>out, submitCmd
    ## Loop over all systematics
    for systName in systAll:
        arg = systAll[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
    for systName in systMC:
        arg = systMC[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))

for d in dataList:
    name = d['name']
    submitCmd  = "create-batch --cfg analyze_data_cfg.py --maxFiles 100"
    submitCmd += " --jobName %s --fileList %s/dataset_%s.txt" % (name, dataDir, name)

    print>>out, submitCmd
    ## Loop over all systematics
    for systName in systAll:
        arg = systAll[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))

os.system("chmod +x %s/submit.sh" % outDir)
