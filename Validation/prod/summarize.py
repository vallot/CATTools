#!/usr/bin/env python

import os
from ROOT import *

out = {}
for d in os.listdir("."):
    if not os.path.isdir(d): continue

    sumW = 0.0
    nTotal = 0
    for fName in os.listdir(d):
        if not fName.endswith(".root"): continue

        f = TFile(d+"/"+fName)
        if f == None: continue
        h = f.Get("gen/hWeight_Norm")
        if h == None: continue
        sumW += h.GetEntries()*h.GetMean()
        nTotal += h.GetEntries()
    if nTotal == 0: out[d] = [0, 0]
    else: out[d] = [nTotal, sumW/nTotal]

maxW = max([len(x) for x in out.keys()])
frm = "%"+str(maxW)+"s  %10d %10.8f"
print "="*(maxW+2+10+2+10)
for name in sorted(out.keys()):
    n, w = out[name]
    if n == 0 or w == 1: continue
    #if n == 0: continue
    print frm % (name, n ,w)

