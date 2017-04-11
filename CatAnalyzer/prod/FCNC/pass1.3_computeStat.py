#!/usr/bin/env python

import os
from glob import glob
from ROOT import *

## Load JSON file and categorize datasets
import json
from pandas import DataFrame
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
dsets = DataFrame.from_dict(json.loads(open("%s/dataset.json" % dataDir).read()))

fsets = DataFrame(columns=('name', 'dsName', 'syst', 'nevt', 'sumW', 'files'))
for path, dirs, files in os.walk('pass2'):
    if len(files) == 0: continue

    ## Find root files and compute basic statistics
    fPaths = []
    nevt, sumW = 0., 0.
    for fName in files:
        if not fName.endswith('.root'): continue
        fPath = path+'/'+fName
        f = TFile(fPath)
        if f == None or f.IsZombie(): continue
        fPaths.append(fPath)

        h = f.Get("agen/hWeight")
        if h == None: h = f.Get("gen/hWeight_Norm")
        if h == None: continue

        nevt += h.GetEntries()
        sumW += h.GetEntries()*h.GetMean()
    if len(fPaths) == 0: continue

    ## Extract stat from the database
    x = path.split('/')
    name, syst = x[1], x[2],
    if len(x) > 3: syst += '/'+x[3]

    dsName = name.split('.')[0]

    fsets.loc[fsets.shape[0]] = [name, dsName, syst, nevt, sumW, fPaths]
    
dsOut = {}
for fsIndex, fset in fsets.iterrows():
    dsName = fset['dsName']
    dset = dsets[dsets.name == dsName]
    if len(dset) != 1:
        print "Multple matching!!!"
        print dset
    dset = dset.iloc[0]

    title, colour, dtype = dset['title'], dset['colour'], dset['type']
    orig_nevt, orig_avgWgt = dset['nevt'], dset['avgWgt']
    lumi, xsec = dset['lumi'], dset['xsec']

    avgWgt = 1
    if fset['nevt'] != 0: avgWgt = fset['sumW']/fset['nevt']
    if "." in fset['name']: title += ":"+" ".join(fset['name'].split('.')[1:])

    newName = '%s/%s' % (fset['syst'], fset['name'])
    dsOut[newName] = {
        'dsName':fset['dsName'],
        'files':fset['files'],
        'nevt':fset['nevt'], 'sumW':fset['sumW'], 'avgWgt':avgWgt,
        'orig_nevt':orig_nevt, 'orig_avgWgt':orig_avgWgt,
        'title':title, 'colour':colour, 'type':dtype,
        'lumi':lumi, 'xsec':xsec,
    }

open("pass2/dataset.json", "w").write(json.dumps(dsOut, sort_keys=True, indent=4))
