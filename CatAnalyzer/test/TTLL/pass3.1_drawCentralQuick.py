#!/usr/bin/env python

import sys, os
import json
sys.argv.append("-b")
from ROOT import *
from CATTools.CatAnalyzer.tdrstyle import *
setTDRStyle()

lumi = 2.11*1000

srcMCs = [
    ["t_bar_t__Jets_LL", 632],
    ["t_bar_t__Jets_Others", 632+1],
    ["Single_top", 800,],
    ["Dibosons", 432,],
    ["Tribosons", 433],
    ["W_Jets", 416],
    ["Z__gamma_rightarrow_ll", 600],
]
for s in srcMCs: s.append(TFile("pass2/central/%s.root" % s[0]))

srcRD_ee = "DoubleEG"
srcRD_mm = "DoubleMuon"
srcRD_em = "MuonEG"

## Load input templates for central result
#sjs = json.loads(open("pass2/samples.json").read())

## Pick the first root file to get full list of plots
plts = []
f = TFile("pass2/central/%s.root" % srcMCs[0][0])
moddir = f.Get("ttll")
for ch in [x.GetName() for x in moddir.GetListOfKeys()]:
    chdir = moddir.GetDirectory(ch)
    if chdir == None: continue

    for step in [x.GetName() for x in chdir.GetListOfKeys()]:
        stepobj = chdir.Get(step)
        if stepobj == None: continue

        if stepobj.IsA().GetName() in ("TH1D", "TH1F"):
            plts.append("ttll/%s/%s" % (ch, step))
        elif stepobj.IsA().InheritsFrom("TDirectory"):
            for plt in [x.GetName() for x in stepobj.GetListOfKeys()]:
                if stepobj.Get(plt) == None: continue
                plts.append("ttll/%s/%s/%s" % (ch, step, plt))

## Prepare output
if not os.path.exists("pass3"): os.mkdir("pass3")

## Save plot list
f = open("pass3/plots.json", "w")
f.write(json.dumps({'plots':plts}, indent=4, sort_keys=True))
f.close()

## Start loop
fout = TFile("pass3/central.root", "recreate")
for plt in plts:
    print "Plotting", plt
    dirName = os.path.dirname(plt)
    pltName = os.path.basename(plt)
    if fout.GetDirectory(dirName) == None: fout.mkdir(dirName)
    fout.cd(dirName)

    c = TCanvas("c_%s" % pltName, pltName, 500, 500)

    ## Add real data histograms
    fRD = None
    if   plt.startswith("ttll/ee/"): fRD = TFile("pass2/central/%s.root" % srcRD_ee)
    elif plt.startswith("ttll/mm/"): fRD = TFile("pass2/central/%s.root" % srcRD_mm)
    elif plt.startswith("ttll/em/"): fRD = TFile("pass2/central/%s.root" % srcRD_em)
    else: print plt

    hRD = fRD.Get(plt).Clone()
    hRD.SetOption("pe")
    hRD.SetMarkerSize(5)

    ## Add MC histograms
    hsMC = THStack("hsMC", "hsMC")
    for finName, color, f in srcMCs:
        h = f.Get(plt)#.Clone()
        h.Scale(lumi)
        h.SetOption("hist")
        h.SetFillColor(color)
        h.SetLineColor(color)
        #h.SetLineStyle(0)
        hsMC.Add(h)
        #del(h)

    ## Draw'em all
    hRD.Draw()
    hsMC.Draw("same")
    hRD.Draw("samep")
    #c.Modified()
    #c.Update()
    fout.cd(dirName)
    c.Write()

    del(hRD)
    del(hsMC)
    del(c)

