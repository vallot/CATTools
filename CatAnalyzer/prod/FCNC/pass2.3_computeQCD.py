#!/usr/bin/env python

import sys, os
import json
#from array import array
from ROOT import *
from math import sqrt

def rootGetDir(srcDir, outPath):
    d = srcDir
    for dName in outPath.split('/'):
        dd = d.GetDirectory(dName)
        if dd == None: dd = d.mkdir(dName)
        d = dd
    return d

## Dataset information
lumi = 36.8*1000
fRD = {
    'el':TFile("pass3/nominal/SingleElectron.root"),
    'mu':TFile("pass3/nominal/SingleMuon.root"),
}
fRD_QCD = {
    'el':TFile("pass3/antiIso/SingleElectron.root"),
    'mu':TFile("pass3/antiIso/SingleMuon.root"),
}

## Load list of datasets, histograms
dataset = json.loads(open("pass3/dataset.json").read())
hists = {
    "el":[x for x in json.loads(open("pass3/hists.json").read()) if x.startswith("el/el")],
    "mu":[x for x in json.loads(open("pass3/hists.json").read()) if x.startswith("mu/mu")],
}

bkgMCs = [
    ["t_bar_t__Jets_rightarrow_l___pm_", 632],
    ["t_bar_t__Jets_Others", 632+3],
    ["SingleTop", 800,],
    ["Dibosons", 432,],
    ["Tribosons", 433],
    ["Z__gamma_rightarrow_ll", 600],
    ["W_Jets", 416],
]

for s in bkgMCs:
    s.append([TFile("pass3/nominal/%s.root" % s[0])])
    if os.path.exists("pass3/antiIso/%s.root" % s[0]): s[-1].append(TFile("pass3/antiIso/%s.root" % s[0]))

## histograms to compute A/B/C/D
dNames = {
    "el":sorted(list(set([os.path.dirname(x) for x in hists["el"]]))),
    "mu":sorted(list(set([os.path.dirname(x) for x in hists["mu"]]))),
}

## Compute scales, 
scales = {}
counts = {}
f = TFile("pass3/nominal/sel_ABCD.root", "recreate")
for ch in ("el", "mu"):
    for dName in dNames[ch]:
        dout = rootGetDir(f, dName)
        ## Make ABCD 
        dout.cd()
        hABCD = TH2D("event_ABCD", "ABCD;isIso;isQCD", 2, 0, 2, 2, 0, 2)

        hName = str(dName+"/event_ABCD")
        hAB, hCD = fRD[ch].Get(hName), fRD_QCD[ch].Get(hName)
        if None in (hAB, hCD): continue

        hABCD.Add(hAB)
        hABCD.Add(hCD)
        for title, colour, files in bkgMCs:
            if len(files) < 2: continue
            hABCD.Add(files[0].Get(hName), -lumi)
            hABCD.Add(files[1].Get(hName), -lumi)
        hABCD.Write()

        ## Set the default scale
        scale = 0

        nA = max(0, hABCD.GetBinContent(2,1))
        nB = max(0, hABCD.GetBinContent(2,2))
        nC = max(0, hABCD.GetBinContent(1,1))
        nD = max(0, hABCD.GetBinContent(1,2))
        if nB > 0 and nD > 0: scale = nB/nD

        scales[dName] = scale
        counts[dName] = scale*(nC+nD)
f.Close()

## Do the subtraction to get the shapes
f = TFile("pass3/nominal/QCD_Data.root", "recreate")
for ch in ("el", "mu"):
    for hName in hists[ch]:
        if hName.endswith("_ABCD"): continue
        dName = str(os.path.dirname(hName))
        if dName not in scales: continue

        scale = scales[dName]
        hName = str(hName)

        hQCD = fRD_QCD[ch].Get(hName)
        if hQCD == None: continue

        rootGetDir(f, dName).cd()
        hQCD = hQCD.Clone()

        for title, colour, files in bkgMCs:
            if len(files) < 2: continue
            hQCD.Add(files[1].Get(hName), -lumi*scale)
        hQCD.Scale(1./lumi)

        hQCD.Write()
f.Close()

## Print out event counts
for dName in sorted(counts.keys()):
    print dName, counts[dName]
