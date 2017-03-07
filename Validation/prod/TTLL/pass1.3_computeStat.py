#!/usr/bin/env python

import os
from ROOT import *

## Load JSON file and categorize datasets
import json
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
js = json.loads(open("%s/dataset.json" % dataDir).read())

dsIn = {}
for x in js: dsIn[x['name']] = x

ds = {}
for fName in sorted(os.listdir('pass1')):
    if not fName.endswith('.root'): continue
    f = TFile('pass1/'+fName)
    if f == None:
        print "!!! root file under %s is invalid" % name
        continue
    d = fName.replace('.root', '')

    if d in dsIn:
        ds[d] = dsIn[d]
    elif '_' in d:
        origName = '_'.join(d.split('_')[:-1])
        ds[d] = dict(dsIn[origName])
        ds[d]['title'] += ":"+d.split('_')[-1]
    if d not in ds:
        print "Cannot find corresponding histogram for", d

    ds[d]['hist'] = fName

    h = f.Get("gen/hWeight_Norm")
    xsec = 1.0
    normFactor = 1.0
    avgWgt = float(ds[d]['avgWgt'])
    nEvent = int(ds[d]['nevt'])
    partNEvent = 0
    if h != None:
        xsec = ds[d]['xsec']
        partAvgWgt = h.GetMean()
        partNEvent = h.GetEntries()
        if partNEvent != nEvent or abs(avgWgt-partAvgWgt) > 1e-5:
            print "!!! Inconsistent stats %s" % d
            print "    new(%10f, %d) orig(%10f, %d)" % (partAvgWgt, partNEvent, avgWgt, nEvent)
        normFactor = h.GetEntries()*avgWgt
    ds[d]['normFactor'] = avgWgt*nEvent
    ds[d]['part_normFactor'] = normFactor
    ds[d]['part_nevt'] = partNEvent

open("pass1/dataset.json", "w").write(json.dumps(ds, sort_keys=True, indent=4))

