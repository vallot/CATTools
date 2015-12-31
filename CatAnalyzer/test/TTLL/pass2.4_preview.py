#!/usr/bin/env python

import sys, os
import json
from array import array
sys.argv.append("-b")
from ROOT import *
from CATTools.CatAnalyzer.tdrstyle import *
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

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
            plts.append({'name':"ttll/%s/%s" % (ch, step)})
        elif stepobj.IsA().InheritsFrom("TDirectory"):
            for plt in [x.GetName() for x in stepobj.GetListOfKeys()]:
                if stepobj.Get(plt) == None: continue
                plts.append({'name':"ttll/%s/%s/%s" % (ch, step, plt)})

## Start loop
fout = TFile("pass2/preview.root", "recreate")
for iplt, pltInfo in enumerate(plts):
    plt = pltInfo['name']
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
    nbinsX = hRD.GetNbinsX()
    hRD.SetOption("pe")
    hRD.SetMarkerSize(5)
    stats = array('d', [0.]*7)
    hRD.GetStats(stats)
    hRD.AddBinContent(nbinsX, hRD.GetBinContent(nbinsX+1))
    hRD.PutStats(stats)

    ## Add MC histograms
    hsMC = THStack("hsMC", "hsMC")
    for finName, color, f in srcMCs:
        h = f.Get(plt)
        h.Scale(lumi)
        h.GetStats(stats)
        h.AddBinContent(nbinsX, h.GetBinContent(nbinsX+1))
        h.PutStats(stats)
        h.SetOption("hist")
        h.SetFillColor(color)
        h.SetLineColor(color)
        #h.SetLineStyle(0)
        hsMC.Add(h)

    ## Draw'em all
    hRD.SetMinimum(0)
    hRD.Draw()
    hsMC.Draw("samehist")
    hRD.Draw("samep")
    fout.cd(dirName)
    c.Write()

    yMax = max([hsMC.GetHistogram().GetBinContent(i) for i in range(1, nbinsX)])
    yMaxR = max([hsMC.GetHistogram().GetBinContent(i) for i in range(nbinsX/2, nbinsX)])
    yMax = max(yMax, max([hRD.GetBinContent(i) for i in range(1, nbinsX)]))
    yMaxR = max(yMaxR, max([hRD.GetBinContent(i) for i in range(nbinsX/2, nbinsX)]))

    plts[iplt]['yMax'] = yMax
    plts[iplt]['yMaxR'] = yMaxR

    del(hRD)
    del(hsMC)
    del(c)

print "A preview root file for the central sample is produced"
print "Run `root -l pass2/preview.root' and browse into each directories to open canvases"
print "You can also use dumpRoot command from the hep-tools"

## Save plot list
f = open("pass2/plots.json", "w")
f.write(json.dumps({'plots':plts}, indent=4, sort_keys=True))
f.close()
