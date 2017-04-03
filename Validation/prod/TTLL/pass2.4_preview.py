#!/usr/bin/env python

import sys, os
import json
from array import array
sys.argv.append("-b")
from math import hypot
from ROOT import *
import imp
printCutflow = imp.load_source("printCutflow", "../submacros/printCutflow.py").printCutflow
#st = imp.load_source("st", "submacros/tdrstyle.py")
gROOT.LoadMacro("../submacros/tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

dataset = json.loads(open("pass2/dataset.json").read())

srcMCs = [
    ["t#bar{t}+Jets#rightarrow l^{+}l^{-}", 632],
    ["t#bar{t}+Jets:Others", 632+3],
    ["t#bar{t}+V", 632+5],
    ["SingleTop", 800],
    ["Dibosons", 432],
    ["Tribosons", 433],
    ["W+Jets", 416],
    ["Z/#gamma#rightarrow ll", 600],
]
for s in srcMCs: s.append(dataset[s[0]]['hist'])
for s in srcMCs: s.append(TFile(s[-1]))

fRD = {
    'ee':TFile("pass2/DoubleEG.root"),
    'mm':TFile("pass2/DoubleMuon.root"),
    'em':TFile("pass2/MuonEG.root"),
}

## Data driven corrections
if not os.path.exists('pass2/scaler_DY.json'): scaleDY = {'scale':{}}
else: scaleDY = json.loads(open("pass2/scaler_DY.json").read())

## Pick the first root file to get full list of plots
plts = []
f = srcMCs[-1][-1]
moddir = f.Get("eventsTTLL")
for ch in [x.GetName() for x in moddir.GetListOfKeys()]:
    chdir = moddir.GetDirectory(ch)
    if chdir == None: continue

    for step in [x.GetName() for x in chdir.GetListOfKeys()]:
        stepobj = chdir.Get(step)
        if stepobj == None: continue

        if stepobj.IsA().GetName() in ("TH1D", "TH1F"):
            plts.append({'name':"eventsTTLL/%s/%s" % (ch, step)})
        elif stepobj.IsA().InheritsFrom("TDirectory"):
            for plt in [x.GetName() for x in stepobj.GetListOfKeys()]:
                if stepobj.Get(plt) == None: continue
                plts.append({'name':"eventsTTLL/%s/%s/%s" % (ch, step, plt)})

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
    mode = plt.split('/')[1]
    if mode not in fRD: continue
    if fRD[mode].Get(plt) == None: continue

    dataName = os.path.basename(fRD[mode].GetName())[:-5]
    lumi = 1000*sum(x['lumi'] for x in dataset[dataName]['subsamples'])

    hRD = fRD[mode].Get(plt).Clone()
    nbinsX = hRD.GetNbinsX()
    hRD.SetOption("pe")
    hRD.SetMarkerSize(0.7)
    hRD.SetMarkerStyle(kFullCircle)
    stats = array('d', [0.]*7)
    hRD.GetStats(stats)
    hRD.AddBinContent(nbinsX, hRD.GetBinContent(nbinsX+1))
    hRD.PutStats(stats)
    hRD.SetTitle("")
    hRD.SetLineColor(kBlack)
    hRD.SetMarkerColor(kBlack)
    
    ## Add MC histograms
    hsMC = THStack("hsMC", "")
    hMC = hRD.Clone()
    hMC.Reset()
    histsMC = []
    for title, color, finName, f in srcMCs:
        h = f.Get(plt)
        h.Scale(lumi)
        if finName == "Z__gamma_rightarrow_ll" and dirName in scaleDY["scale"]:
            h.Scale(scaleDY["scale"][dirName])
        hMC.Add(h)
        h.GetStats(stats)
        h.AddBinContent(nbinsX, h.GetBinContent(nbinsX+1))
        h.PutStats(stats)
        h.SetOption("hist")
        h.SetFillColor(color)
        h.SetLineColor(color)
        #h.SetLineStyle(0)
        h.GetXaxis().SetLabelOffset(999)
        h.GetXaxis().SetLabelSize(0)
        hsMC.Add(h)
        h.SetTitle(title)
        histsMC.append(h)
    hRatio = hRD.Clone()
    hRatio.Reset()
    hRatio.SetTitle(";%s;Data/MC" % hRD.GetXaxis().GetTitle())
    grpRatio = TGraphErrors()
    rMin, rMax = 0.5, 1.5
    for b in range(nbinsX):
        yRD, yMC = hRD.GetBinContent(b+1), hMC.GetBinContent(b+1)
        eRD, eMC = hRD.GetBinError(b+1), hMC.GetBinError(b+1)
        r, e = 1e9, 1e9
        if yMC > 0:
            r = yRD/yMC
            #rMax = max(r, rMax)
        #if yMC > 0 and yRD > 0: e = r*hypot(eRD/yRD, eMC/yMC)
        if yRD > 0: e = r*abs(eRD/yRD)
        if r == 1e9 or e == 1e9: continue

        x = hRD.GetXaxis().GetBinCenter(b+1)
        w = hRD.GetXaxis().GetBinWidth(b+1)
        n = grpRatio.GetN()
        grpRatio.SetPoint(n, x, r)
        grpRatio.SetPointError(n, w/2, e)
    #if rMax > 2: rMax = 3
    hRatio.SetStats(False)
    hRatio.SetMinimum(rMin)
    hRatio.SetMaximum(rMax)

    hRD.GetXaxis().SetLabelOffset(999)
    hRD.GetXaxis().SetLabelSize(0)

    ## Draw'em all
    plotDim = (400, 300, 100) # width, main height, ratio height
    margin = (2*40, 20, 60, 60) # left, right, bottom, top
    padH = (plotDim[1] + margin[3], plotDim[2] + margin[2])
    canH = padH[0] + padH[1]
    canW = plotDim[0] + margin[0] + margin[1]

    leg = TLegend(0.6, 0.6, 1.-1.*margin[1]/canW-.03, 1-1.*margin[3]/padH[0]-.05)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for h in reversed(histsMC):
        leg.AddEntry(h, h.GetTitle(), "f")
    leg.AddEntry(hRD, "Data", "lp")

    hRatio.GetXaxis().SetTitleSize(0.1)
    hRatio.GetXaxis().SetTitleOffset(0.75)
    hRatio.GetXaxis().SetLabelSize(0.08)

    hRatio.GetYaxis().SetTitleSize(0.1)
    hRatio.GetYaxis().SetTitleOffset(0.75)
    hRatio.GetYaxis().SetLabelSize(0.1)
    hRatio.GetYaxis().SetNdivisions(505)

    hRD.SetStats(False)
    hRD.GetXaxis().SetTitle("")

    hRD.GetYaxis().SetTitleSize(0.045)
    hRD.GetYaxis().SetTitleOffset(1.55)
    hRD.GetYaxis().SetLabelSize(0.045)

    c = TCanvas("c_%s" % pltName, pltName, canH, canW)
    c.Divide(1,2)

    pad2 = c.cd(2)
    pad2.SetPad(0, 0, 1, 1.0*padH[1]/canH)
    pad2.SetMargin(1.*margin[0]/canW, 1.*margin[1]/canW, 1.*margin[2]/padH[1], 0)
    hRatio.Draw()
    grpRatio.Draw("P")
    grpRatio.SetMarkerStyle(kFullCircle)
    grpRatio.SetMarkerSize(0.7)
    pad2.RedrawAxis()

    pad1 = c.cd(1)
    pad1.SetPad(0, 1.0*padH[1]/canH, 1, 1)
    pad1.SetMargin(1.*margin[0]/canW, 1.*margin[1]/canW, 0.05, 1.*margin[3]/padH[0])

    pad1.SetLogy()
    hRD.SetMaximum(hRD.GetMaximum()*100)
    hRD.SetMinimum(1e-3)
    hRD.Draw("")
    hsMC.Draw("samehist")
    hRD.Draw("same,e")
    leg.Draw()
    pad1.RedrawAxis()

    fout.cd(dirName)
    c.Write()

    yMax = max([hsMC.GetHistogram().GetBinContent(i) for i in range(1, nbinsX)])
    yMaxR = max([hsMC.GetHistogram().GetBinContent(i) for i in range(nbinsX/2, nbinsX)])
    yMax = max(yMax, max([hRD.GetBinContent(i) for i in range(1, nbinsX)]))
    yMaxR = max(yMaxR, max([hRD.GetBinContent(i) for i in range(nbinsX/2, nbinsX)]))

    plts[iplt]['yMax'] = yMax
    plts[iplt]['yMaxR'] = yMaxR

    if not os.path.exists("preview/%s" % dirName):
        os.makedirs("preview/%s" % dirName)
    c.Print("preview/%s/%s.png" % (dirName, c.GetName()))
    #if grpRatio.GetN() > 0: c.Print("preview/%s/%s.C" % (dirName, c.GetName()))

    for h in (hRD, hMC, hsMC, hRatio, grpRatio, c): del(h)

## Start to print cut flow
cutflow = {
    "count":{"ee":{}, "mm":{}, "em":{}},
    "error":{"ee":{}, "mm":{}, "em":{}},
    "nstep":0,
    "step":None,
}
nstep = 0
for mode in cutflow["count"].keys():
    h = fRD[mode].Get("eventsTTLL/%s/cutstep" % mode)
    nstep = h.GetNbinsX()
    cutflow["count"][mode]["Data"] = [h.GetBinContent(i) for i in range(1, nstep+1)]
    cutflow["error"][mode]["Data"] = [h.GetBinError(i) for i in range(1, nstep+1)]
    if cutflow["step"] == None:
        cutflow["step"] = [h.GetXaxis().GetBinLabel(i) for i in range(1, nstep+1)]

    for title, color, finName, f in srcMCs:
        h = f.Get("eventsTTLL/%s/cutstep" % mode)
        cutflow["count"][mode][title] = [h.GetBinContent(i) for i in range(1, nstep+1)]
        cutflow["error"][mode][title] = [h.GetBinError(i) for i in range(1, nstep+1)]
cutflow["nstep"] = nstep
printCutflow(cutflow)

## Save cut flow
f = open("pass2/cutflow.json", "w")
f.write(json.dumps(cutflow, indent=4, sort_keys=True))
f.close()

## Save plot list
f = open("pass2/plots.json", "w")
f.write(json.dumps({'plots':plts}, indent=4, sort_keys=True))
f.close()

print "A preview root file is produced"
print "Run `root -l pass2/preview.root' and browse into each directories to open canvases"
print "You can also use dumpRoot command from the hep-tools, dumpRoot pass2/preview.root"
print "Cutflow table is saved in JSON format, pass2/cutflow.json"
