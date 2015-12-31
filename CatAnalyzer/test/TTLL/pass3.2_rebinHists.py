#!/usr/bin/env python

import sys, os
import json
from array import array
from numpy import linspace
sys.argv.append("-b")
from ROOT import *

def arr(min, max, n = None):
    if n == None: n = max-min
    return array('d', list(linspace(min, max, n+1)))

plotInfos = json.loads(open("pass3/plots.json").read())
binInfos = {
    "vertex_n":arr(0, 51),
    #"met_pt":array('d', [0, ]),
    #"met_phi":array('d', []),
    "jets_n":arr(0,6),
    "bjets_n":arr(0,6),
    "lepton1_pt":arr(0,201),
    "lepton2_pt":arr(0,201),
    "lepton1_eta":arr(-2.5, 2.5, 50),
    "lepton2_eta":arr(-2.5, 2.5, 50),
    "z_m":(array('d', [0, 20, 30, 40, 50, 60]+list(linspace(91-15, 91+15, 5))+[120, 150, 200])),
}

for fiPath, dirs, files in os.walk("pass2"):
    rootFiles = [x for x in files if x.endswith('.root')]
    if len(rootFiles) == 0: continue

    foPath = fiPath.replace("pass2", "pass3/hists")
    if not os.path.exists(foPath): os.makedirs(foPath)

    for fName in rootFiles:
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
            ho = hi.Clone()
            if hName in binInfos: 
                bins = binInfos[hName]
                ho = ho.Rebin(len(bins)-1, "", bins)
            ho.Write()
