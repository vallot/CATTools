#!/usr/bin/env python

import sys, os
outDir = "pass1_makeHist"
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

## Write script to run create-batch
out = open("%s/submit_sig.sh" % outDir, "w")
for d in sigList:
    name = d['name']
    title = d['title']
    submitCmd  = "create-batch --cfg analyze_sig_cfg.py --maxFiles 25"
    submitCmd += " --fileList %s/dataset_%s.txt" % (dataDir, name)

    print>>out, (submitCmd + " --jobName %s/central" % name)

    ## Scale up/down systematic undertainty from LHE weight
    ## This uncertainty have to be combined with envelope
    ## Let us assume index1-10 are for the scale variations (muF & muR)
    if '_aMC' in name or '_powheg' in name:
        print>>out, "## Scale variations in aMC@NLO sample"
        for i in range(1,9): # total 8 scale variations, 3 muF x 3 muR and one for central weight
            arg =  'src="genWeight:pdfWeights" genWeight.index=%d' % i
            print>>out, (submitCmd + (" --jobName %s/gen_scale/%d --args '%s'" % (name, i, arg)))

        if '_scale' not in name:
            if '_aMC' in name: weightSize = 110
            elif '_powheg' in name: weightSize = 248
            ## NOTE: there is weight vector, but we don't do it for LO generator here.
            ##elif '_madgraph' in name: weightSize = 445
            else: weightSize = 0
            for i in range(9, weightSize+1):
                arg = 'src="genWeight:pdfWeights" genWeight.index=%d' % i
                print>>out, (submitCmd + (" --jobName %s/gen_PDF/%d --args '%s'" % (name, i, arg)))

    ## Loop over all systematics
    for systName in systAll:
        arg = systAll[systName]
        print>>out, (submitCmd + (" --jobName %s/%s --args '%s'" % (name, systName, arg)))
    for systName in systMC:
        arg = systMC[systName]
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
