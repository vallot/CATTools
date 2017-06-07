#!/usr/bin/env python

import os
from ROOT import *

## Load JSON file and categorize datasets
import json
from pandas import DataFrame
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
dsets = DataFrame.from_dict(json.loads(open("%s/dataset.json" % dataDir).read()))

fsets = DataFrame(columns=('hist', 'name', 'nevt', 'avgWgt', 'options'))
for fName in os.listdir('pass1'):
    if not fName.endswith(".root"): continue
    f = TFile('pass1/'+fName)
    if f == None or f.IsZombie(): continue

    tokens = fName[:-5].split('.')
    name = tokens.pop(0)

    h = f.Get("gen/hWeight_Norm")
    nevt, avgWgt = 0, 1
    if h != None:
        nevt = h.GetEntries()
        avgWgt = h.GetMean()

    fsets.loc[fsets.shape[0]] = [fName, name, nevt, avgWgt, tokens[:]]

dsOut = {}
for dsIndex, dsRow in dsets.iterrows():
    name = dsRow['name']
    fsRows = fsets[fsets['name'].str.match('^'+name+'$')]

    for fsIndex, fsRow in fsRows.iterrows():
        fName = fsRow['hist']
        newName = fName.split('/')[-1][:-5]

        nevt, avgWgt = 0.0, 1.0
        nevt_orig, avgWgt_orig = dsRow.nevt, dsRow.avgWgt
        normFactor = 1.0

        if dsRow.type != 'Data':
            nevt += fsRow['nevt']
            avgWgt += fsRow['avgWgt']*fsRow['nevt']
        if nevt > 0: avgWgt /= nevt
        normFactor = nevt*avgWgt

        dsOut[newName] = {}
        dsOut[newName].update(dict(dsRow))
        dsOut[newName].update({'hist':'pass1/'+fName, 'nevt':nevt, 'avgWgt':avgWgt, 'normFactor':normFactor})
        if len(fsRow['options']) > 0:
            dsOut[newName]['title'] += ':'+('+'.join(fsRow['options']))

        print newName, nevt, avgWgt

open("pass1/dataset.json", "w").write(json.dumps(dsOut, sort_keys=True, indent=4))
