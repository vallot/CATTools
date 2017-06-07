#!/usr/bin/env python

import sys, os
if not os.path.exists("pass3"): os.mkdir("pass3")

from math import sqrt
import json
from ROOT import *
def findHists(d, hists):
    if d == None: return
    for name in [key.GetName() for key in d.GetListOfKeys()]:
        obj = d.Get(name)
        if obj == None: continue
        if obj.IsA().InheritsFrom("TDirectory"):
            findHists(obj, hists)
        else:
            hPath = (d.GetPath().split(':')[-1])+"/"+obj.GetName()
            hists.add(hPath[1:])

checkedTypes = set()
ds = {}
dsSyst = {}
dsIn = json.loads(open("pass2/dataset.json").read())
hists = set()
for name in dsIn:
    x = dsIn[name]
    syst = os.path.dirname(name)
    title = "%s/%s" % (syst, x["title"])
    stype = x["type"]
    safeTitle = x["title"]
    for i in " ;:!@#$%^&*()-+=/<>?[]{}|": safeTitle = safeTitle.replace(i,'_')
    if title not in ds:
        ds[title] = {
            'colour':x['colour'],
            'hist':'pass3/%s/%s.root' % (syst, safeTitle), ## Path to the merged histogram
            'subsamples':[], ## List of input samples
        }
    scale = 1
    if stype != 'Data' and x['sumW'] != 0: scale = x['xsec']/x['sumW']
    if scale == 0: scale = 1

    ds[title]['subsamples'].append({
        'type':stype,
        'avgWgt':x['avgWgt'], 'nevt':x['nevt'],
        'xsec':x['xsec'], 'lumi':x['lumi'], 'sumW':x['sumW'], 'scale':scale,
        'files':x['files'],
    })

    sstype = stype.split('/')[0]
    if sstype not in checkedTypes or len(hists) == 0:
        print "Reading histogram contents of sample type", sstype
        hset = set()
        f = TFile(x['files'][0])
        findHists(f, hset)
        f.Close()
        hists = hists.union(hset)
        checkedTypes.add(sstype)

hists = list(hists)
hists.sort()

open("pass3/dataset.json", "w").write(json.dumps(ds, sort_keys=True, indent=4))
open("pass3/hists.json", "w").write(json.dumps(hists, sort_keys=True, indent=4))

