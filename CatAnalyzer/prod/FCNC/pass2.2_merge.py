#!/usr/bin/env python
import sys, os
import json
from ROOT import *

def rootmkdirs(f, path):
    d = f.GetDirectory(path)
    if d != None: return d

    d = f
    for t in path.split('/'):
        td = d.GetDirectory(t)
        if td != None: d = td
        else: d = d.mkdir(t)
    return d

def mergeHist(dout, din, scale, allHists):
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
        mergeHist(dNext, din.GetDirectory(dName), scale, allHists)

    for hName in hNames:
        hPath = dout.GetPath()+"/"+hName
        if hPath not in allHists:
            dout.cd()
            h = din.Get(hName).Clone()
            h.Reset()
            allHists[hPath] = h

        h = allHists[hPath]
        h.Add(din.Get(hName), scale)

ds = json.loads(open("pass3/dataset.json").read())
hs = json.loads(open("pass3/hists.json").read())

for d in ds:
    print "Processing dataset", d
    foutName = ds[d]['hist']
    outPath = os.path.dirname(foutName)

    if not os.path.isdir(outPath): os.makedirs(outPath)

    sses = []
    for ss in ds[d]['subsamples']:
        scale = ss['scale']
        sses.append( ([TFile(fName) for fName in ss['files']], scale) )

    fout = TFile(foutName, "RECREATE")
    for hName in hs:
        hName = str(hName)
        ## Check the histogram exists in the source root file
        if sses[0][0][0].Get(hName) == None: continue

        ## Prepare output directory
        p = '/'.join(hName.split('/')[:-1])
        dout = rootmkdirs(fout, p)
        dout.cd()

        ## Prepare output histogram
        h = sses[0][0][0].Get(hName).Clone()
        h.SetDirectory(dout)
        h.Reset()
        #h.Sumw2()
        ## Merge histograms
        for files, scale in sses:
            for fin in files:
                hin = fin.Get(hName)
                if not hin.IsA().InheritsFrom("TH1"): continue
                h.Add(hin, scale)
                hin.Delete()
        dout.cd()
        h.Write()

    for ss in sses:
        for f in ss[0]: f.Close()
    sses = None
    fout.Close()
