#!/usr/bin/env python
import sys, os
import json
from ROOT import *

def mergeHist(dout, din, normFactor, allHists):
    if din == None or dout == None: return allHists

    dNames, hNames = [], []
    for name in [key.GetName() for key in din.GetListOfKeys()]:
        obj = din.Get(name)
        if obj == None: continue

        if obj.IsA().InheritsFrom("TDirectory"): dNames.append(name)
        elif obj.IsA().InheritsFrom("TH1"): hNames.append(name)

    for dName in dNames:
        dNext = dout.GetDirectory(dName)
        if dNext == None: dNext = dout.mkdir(dName)
        mergeHist(dNext, din.GetDirectory(dName), normFactor, allHists)

    for hName in hNames:
        hPath = dout.GetPath()+"/"+hName
        if hPath not in allHists:
            dout.cd()
            h = din.Get(hName).Clone()
            h.Scale(1./normFactor)
            allHists[hPath] = h
        else:
            h = allHists[hPath]
            h.Add(din.Get(hName), 1./normFactor)

ds = json.loads(open("pass2/dataset.json").read())
hists = json.loads(open("pass2/hists.json").read())
for d in ds:
    print "Processing dataset", d
    foutName = ds[d]['hist']
    outPath = os.path.dirname(foutName)
    if not os.path.isdir(outPath): os.mkdir(outPath)
    fout = TFile(foutName, "RECREATE")

    ## Merge histograms
    allHists = {}
    for sample in ds[d]['samples']:
        fInName = sample['hist']
        if not os.path.exists(fInName): continue
        fin = TFile(fInName)
        if fin == None: continue
        print sample

        normFactor = sample['normFactor']
        mergeHist(fout, fin, normFactor, allHists)
        fin.Close()

        for path in allHists:
            fout.GetDirectory('/'.join(path.split('/')[:-1])).cd()
            allHists[path].Write()


