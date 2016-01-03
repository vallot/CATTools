#!/usr/bin/env python

import sys, os
import json
from array import array
sys.argv.append("-b")
from math import hypot
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

fRDs = {
    "ee":TFile("pass2/central/DoubleEG.root"),
    "mm":TFile("pass2/central/DoubleMuon.root"),
    "em":TFile("pass2/central/MuonEG.root"),
}

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

    ## Add real data histograms
    fRD = fRDs[plt[5:7]]

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
    hMC = hRD.Clone()
    hMC.Reset()
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
        hMC.Add(h)
    grpRatio = TGraphErrors()
    grpRatio.SetTitle("ratio;%s;Data/MC" % hRD.GetXaxis().GetTitle())
    rMax = 3
    for b in range(1, nbinsX+1):
        yRD, yMC = hRD.GetBinContent(b), hMC.GetBinContent(b)
        eRD, eMC = hRD.GetBinError(b), hMC.GetBinError(b)
        r, e = 1e9, 1e9
        if yMC > 0:
            r = yRD/yMC
        #    rMax = max(r, rMax)
        if yMC > 0 and yRD > 0: e = r*hypot(eRD/yRD, eMC/yMC)

        x = hRD.GetXaxis().GetBinCenter(b)
        w = hRD.GetXaxis().GetBinWidth(b)
        grpRatio.SetPoint(b, x, r)
        grpRatio.SetPointError(b, w/2, e)
    grpRatio.SetMinimum(0)
    grpRatio.SetMaximum(rMax)

    ## Draw'em all
    c = TCanvas("c_%s" % pltName, pltName, 500, 700)
    c.Divide(1,2)

    pad2 = c.cd(2)
    pad2.SetPad(0, 0, 1, 2./7)
    pad2.SetTopMargin(0)
    grpRatio.GetYaxis().SetLabelSize(grpRatio.GetYaxis().GetLabelSize()*5./2)
    grpRatio.GetYaxis().SetTitleSize(grpRatio.GetYaxis().GetTitleSize()*5./2)
    grpRatio.GetXaxis().SetLabelSize(grpRatio.GetXaxis().GetLabelSize()*5./2)
    grpRatio.GetXaxis().SetTitleSize(grpRatio.GetXaxis().GetTitleSize()*5./2)
    grpRatio.Draw("AP")

    pad1 = c.cd(1)
    pad1.SetPad(0, 2./7, 1, 1)
    pad1.SetBottomMargin(0)

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

## Start to print cut flow
cutflow = {"ee":{}, "mm":{}, "em":{}}
nstep = 0
fws = [0,]
fws[0] = max([len(x[0]) for x in srcMCs]+[4,])
for mode in cutflow.keys():
    h = fRDs[mode].Get("ttll/%s/cutstep" % mode)
    nstep = h.GetNbinsX()
    cutflow[mode]["Data"] = [h.GetBinContent(i) for i in range(1, nstep+1)]
    if len(fws) == 1: fws.extend([len(h.GetXaxis().GetBinLabel(i)) for i in range(1, nstep+1)])

    for finName, color, f in srcMCs:
        h = f.Get("ttll/%s/cutstep" % mode)
        cutflow[mode][finName] = [h.GetBinContent(i) for i in range(1, nstep+1)]
fwtot = sum(fws)+12+len(fws)
for mode in cutflow.keys():
    print "="*((fwtot-12)/2), "Cutflow for", mode, "="*((fwtot-12)/2)
    print " "*fws[0], "|",
    print " | ".join([h.GetXaxis().GetBinLabel(i+1) for i in range(nstep)])
    tfmt = "%"+str(fws[0])+"s |"
    cutflow_bkg = [0.]*nstep
    for x in cutflow[mode]:
        if x == "Data" or 't_bar_t' in x: continue
        print tfmt % x,
        print " | ".join([("%"+str(fws[i+1])+".2f") % cutflow[mode][x][i] for i in range(nstep)])
        for i in range(nstep): cutflow_bkg[i] += cutflow[mode][x][i]
    print "-"*fwtot
    cutflow_sig = [0.]*nstep
    for x in cutflow[mode]:
        if 't_bar_t' not in x: continue
        print tfmt % x,
        print " | ".join([("%"+str(fws[i+1])+".2f") % cutflow[mode][x][i] for i in range(nstep)])
        for i in range(nstep): cutflow_sig[i] += cutflow[mode][x][i]
    print "-"*fwtot
    print tfmt % "All Signal",
    print " | ".join([("%"+str(fws[i+1])+".2f") % cutflow_sig[i] for i in range(nstep)])
    print tfmt % "All Bkg",
    print " | ".join([("%"+str(fws[i+1])+".2f") % cutflow_bkg[i] for i in range(nstep)])
    print tfmt % "All MC",
    print " | ".join([("%"+str(fws[i+1])+".2f") % (cutflow_sig[i] + cutflow_bkg[i]) for i in range(nstep)])
    print "-"*fwtot
    print tfmt % "Data",
    print " | ".join([("%"+str(fws[i+1]-3)+"d   ") % cutflow[mode]["Data"][i] for i in range(nstep)])
    print "="*fwtot
    print

print "A preview root file for the central sample is produced"
print "Run `root -l pass2/preview.root' and browse into each directories to open canvases"
print "You can also use dumpRoot command from the hep-tools, dumpRoot pass2/preview.root"

## Save plot list
f = open("pass2/plots.json", "w")
f.write(json.dumps({'plots':plts}, indent=4, sort_keys=True))
f.close()
