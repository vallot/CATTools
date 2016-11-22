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
hs = json.loads(open("pass2/hists.json").read())

for d in ds:
    print "Processing dataset", d
    foutName = ds[d]['hist']
    outPath = os.path.dirname(foutName)
    if not os.path.isdir(outPath): os.mkdir(outPath)

    sses = []
    for ss in ds[d]['subsamples']:
        fName = ss['hist']
        xsec, normFactor = ss['xsec'], ss['normFactor']
        scale = 1
        if ss['type'] != 'Data' and normFactor != 0: scale = xsec/normFactor
        sses.append( (TFile(fName), scale) )

    fout = TFile(foutName, "RECREATE")
    for hName in hs:
        hName = str(hName)
        ## Check the histogram exists in the source root file
        if sses[0][0].Get(hName) == None: continue

        ## Prepare output directory
        p = '/'.join(hName.split('/')[:-1])
        dout = rootmkdirs(fout, p)
        dout.cd()

        ## Prepare output histogram
        h = sses[0][0].Get(hName).Clone()
        h.SetDirectory(dout)
        h.Reset()
        h.Sumw2()
        ## Merge histograms
        for ss in sses:
            hin = ss[0].Get(hName)
            h.Add(hin, ss[1])
            hin.Delete()
        dout.cd()
        h.Write()

    for ss in sses: ss[0].Close()
    sses = None
    fout.Close()

"""
    ## Merge histograms
    allHists = {}
    for sample in ds[d]['samples']:
        fInName = sample['hist']
        if not os.path.exists(fInName): continue
        fin = TFile(fInName)
        if fin == None: continue

        normFactor = sample['normFactor']
        mergeHist(fout, fin, normFactor, allHists)
        fin.Close()

        for path in allHists:
            fout.GetDirectory('/'.join(path.split('/')[:-1])).cd()
            allHists[path].Write()

"""
