#!/usr/bin/env python

import os
from ROOT import *

## Load JSON file and categorize datasets
import json
import pandas as pd
pd.options.display.max_rows = 999
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
dsets = pd.DataFrame.from_dict(json.loads(open("%s/dataset.json" % dataDir).read()))

fsets = pd.DataFrame(columns=('hist', 'name', 'nevt', 'avgWgt'))
for fName in os.listdir('pass1'):
    if not fName.endswith(".root"): continue
    f = TFile('pass1/'+fName)
    if f == None or f.IsZombie(): continue

    name = fName.split('.')[0]

    h = f.Get("gen/hWeight_Norm")
    nevt, avgWgt = 0, 1
    if h != None:
        nevt = h.GetEntries()
        avgWgt = h.GetMean()

    fsets.loc[fsets.shape[0]] = [fName, name, nevt, avgWgt]

dsets['nevt_new'] = [0.0]*dsets.shape[0]
dsets['avgWgt_new'] = [1.0]*dsets.shape[0]
dsOut = {}
for dsIndex, dsRow in dsets.iterrows():
    name = dsRow['name']
    fsRows = fsets[fsets['name'].str.match('^'+name+'$')]
    if fsRows.shape[0] == 0: continue

    nevt, avgWgt = 0.0, 1.0
    nevt_orig, avgWgt_orig = dsRow.nevt, dsRow.avgWgt
    normFactor = 1.0
    if dsRow.type != 'Data':
        for fsIndex, fsRow in fsRows.iterrows():
            nevt += fsRow['nevt']
            avgWgt += fsRow['avgWgt']*fsRow['nevt']
        if nevt > 0: avgWgt /= nevt
        normFactor = nevt*avgWgt

    dsOut[name] = {}
    dsOut[name].update(dict(dsRow))
    dsOut[name].update({'nevt':nevt, 'avgWgt':avgWgt, 'normFactor':normFactor})

    dsets.loc[dsIndex, 'nevt_new'] = nevt
    dsets.loc[dsIndex, 'avgWgt_new'] = avgWgt

subds = dsets[['name', 'nevt', 'nevt_new', 'avgWgt', 'avgWgt_new']]
print subds[(subds['nevt'] != subds['nevt_new']) | (abs(subds['avgWgt']-subds['avgWgt_new'])>1e-5)]

open("pass1/dataset.json", "w").write(json.dumps(dsOut, sort_keys=True, indent=4))
