#!/usr/bin/env python

import sys, os
import json
from array import array
from numpy import linspace
from multiprocessing import Pool, cpu_count
sys.argv.append("-b")
from ROOT import *

def arr(xmin, xmax, n = None):
    if n == None: n = xmax-xmin
    return array('d', list(linspace(xmin, xmax, n+1)))

plotInfos = json.loads(open("pass2/plots.json").read())
binInfos = {
    "vertex_n":arr(0, 50),
    "met_pt":arr(0,200,20),
    "met_phi":arr(-3.15,3.15,20),
    "jets_n":arr(0,10),
    "bjets_n":arr(0,6),
    "lepton1_pt":arr(0,200),
    "lepton2_pt":arr(0,200),
    "lepton1_eta":arr(-2.5, 2.5, 50),
    "lepton2_eta":arr(-2.5, 2.5, 50),
    "jet1_pt":arr(0,200,20),
    "jet2_pt":arr(0,200,20),
    "jet3_pt":arr(0,200,20),
    "jet4_pt":arr(0,200,20),
    "jet1_eta":arr(-2.5,2.5,20),
    "jet2_eta":arr(-2.5,2.5,20),
    "jet3_eta":arr(-2.5,2.5,20),
    "jet4_eta":arr(-2.5,2.5,20),
    "z_m":arr(20, 320, 50),
}

def process(fiPath, foPath, fName):
    print "Processing", foPath+"/"+fName
    fi = TFile(fiPath+"/"+fName)
    fo = TFile(foPath+"/"+fName, "RECREATE")
    for plotInfo in plotInfos["plots"]:
        hPath = plotInfo["name"]
        hName = os.path.basename(hPath)
        hPath = os.path.dirname(hPath)

        d = fo.GetDirectory(hPath)
        if d == None: fo.mkdir(hPath)
        d = fo.GetDirectory(hPath)
        d.cd()

        hi = fi.Get(hPath+"/"+hName)
        if hi == None: continue
        ho = hi.Clone()
        if hName in binInfos:
            bins = binInfos[hName]
            ho = ho.Rebin(len(bins)-1, "", bins)
        ho.Write()

if __name__ == '__main__':
    p = Pool(cpu_count())
    for fiPath, dirs, files in os.walk("pass2"):
        rootFiles = [x for x in files if x.endswith('.root')]
        if len(rootFiles) == 0: continue

        foPath = fiPath.replace("pass2", "pass3/hists")
        if not os.path.exists(foPath): os.makedirs(foPath)

        for fName in rootFiles:
            p.apply_async(process, [fiPath ,foPath, fName])
    p.close()
    p.join()
