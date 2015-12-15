#!/usr/bin/env python

import sys, os
outDir = "step1_makeHist"
if os.path.exists(outDir):
    print "Directory already exists. remove or rename it."
    sys.exit()
os.mkdir(outDir)
os.system("ln -s ../analyze_sig_cfg.py %s/analyze_sig_cfg.py" % outDir)
os.system("ln -s ../analyze_bkg_cfg.py %s/analyze_bkg_cfg.py" % outDir)
os.system("ln -s ../analyze_data_cfg.py %s/analyze_data_cfg.py" % outDir)

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
    if datasetName.endswith("/MINIAOD"):
        if datasetName.startswith("/Jet"): continue
        dataList.append(j)
    elif title.startswith('t#bar{t}'): sigList.append(j)
    else:
        if datasetName.startswith("/QCD"): continue
        bkgList.append(j)

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
out = open("%s/submit_sig.sh" % outDir, "w")
for d in sigList:
    name = d['name']
    title = d['title']
    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 25"
    submitCmd += " --fileList %s/dataset_%s.txt" % (dataDir, name)

    ## Special care for systematic study samples
    if '_scale' in name:
        continue

    print>>out, (submitCmd + " --jobName %s/central" % name)
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
out.close()

out = open("%s/submit_bkg.sh" % outDir, "w")
for d in bkgList:
    name = d['name']
    submitCmd  = "create-batch --cfg analyze_bkg_cfg.py --maxFiles 25"
    submitCmd += " --jobName %s --fileList %s/dataset_%s.txt" % (name, dataDir, name)

    print>>out, (submitCmd + " --jobName %s/central" % name)
    ## Loop over all systematics
    for systName in systAll:
        arg = systAll[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
    for systName in systMC:
        arg = systMC[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
out.close()

out = open("%s/submit_data.sh" % outDir, "w")
for d in dataList:
    name = d['name']
    submitCmd  = "create-batch --cfg analyze_data_cfg.py --maxFiles 25"
    submitCmd += " --jobName %s --fileList %s/dataset_%s.txt" % (name, dataDir, name)

    print>>out, (submitCmd + " --jobName %s/central" % name)
    ## Loop over all systematics
    for systName in systAll:
        arg = systAll[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
out.close()

os.system("chmod +x %s/submit_sig.sh" % outDir)
os.system("chmod +x %s/submit_bkg.sh" % outDir)
os.system("chmod +x %s/submit_data.sh" % outDir)
