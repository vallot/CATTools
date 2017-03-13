#!/usr/bin/env python

import os
from ROOT import *

## Load JSON file and categorize datasets
import json
from pandas import DataFrame
dataDir = "%s/src/CATTools/CatAnalyzer/data/dataset" % os.environ["CMSSW_BASE"]
dsets = DataFrame.from_dict(json.loads(open("%s/dataset.json" % dataDir).read()))

dsets['hist'] = [None]*len(dsets.name)
dsets['part_nevt'] = [0.0]*len(dsets.name)
dsets['part_avgWgt'] = [1.0]*len(dsets.name)
dsets['normFactor'] = [1.0]*len(dsets.name)

dsOut = {}
for index, name in enumerate(dsets.name):
    fName = 'pass1/%s.root' % name
    if not os.path.exists(fName):
        print '!!!', name, 'cannot be found'
        continue
    f = TFile(fName)
    if f == None or f.IsZombie():
        print '!!!', name, 'cannot opened'
        continue
    dsets.loc[index, 'hist'] = fName

    dsets.loc[index, 'part_nevt'] = dsets.loc[index, 'nevt']
    h = f.Get("gen/hWeight_Norm")
    if h != None:
        dsets.loc[index, 'part_nevt'] = h.GetEntries()
        dsets.loc[index, 'part_avgWgt'] = h.GetMean()
        dsets.loc[index, 'normFactor'] = h.GetEntries()*h.GetMean()

    dsOut[name] = dict(dsets.ix[index]).copy()

subds = dsets[['name', 'nevt', 'part_nevt', 'avgWgt', 'part_avgWgt']]
print subds[(subds.nevt != subds.part_nevt)]


open("pass1/dataset.json", "w").write(json.dumps(dsOut, sort_keys=True, indent=4))
