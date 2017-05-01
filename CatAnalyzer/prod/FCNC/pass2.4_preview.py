#!/usr/bin/env python

import sys, os
import json
from array import array
sys.argv.append("-b")
from math import hypot
from ROOT import *
import imp
printCutflow = imp.load_source("printCutflow", "submacros/printCutflow.py").printCutflow
#st = imp.load_source("st", "submacros/tdrstyle.py")
gROOT.LoadMacro("submacros/tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

variation = "nominal"
#variation = "antiIso"
lumi = 35.87*1000

bkgMCs = [
    ["t_bar_t__Jets_rightarrow_l___pm_", "t#bar{t}+Jets#rightarrow l^{#pm}", 632],
    ["t_bar_t__Jets_Others", "t#bar{t}+Jets Others", 632+3],
    ["SingleTop", "Single top", 800,],
    ["Dibosons", "Dibosons", 432,],
    ["Tribosons", "Tribosons", 433],
    ["Z__gamma_rightarrow_ll", "Z/#gamma#rightarrow ll", 600],
    ["W_Jets", "W+Jets", 416],
    ["QCD_Data", "QCD", kGray],
]
sigMCs = [
    ["tH_rightarrow_t_bar_t_c", "tcH", kBlue],
    ["tH_rightarrow_t_bar_t_u", "tuH", kBlue],
]
for s in bkgMCs+sigMCs:
    if os.path.exists("pass3/%s/%s.root" % (variation, s[0])):
        s.append(TFile("pass3/%s/%s.root" % (variation, s[0])))
    else:
        s.append(TFile("pass3/%s/%s.root" % ("nominal", s[0])))

fRD = {
    'el':TFile("pass3/%s/SingleElectron.root" % variation),
    'mu':TFile("pass3/%s/SingleMuon.root" % variation),
}

## Data driven corrections
dataset = json.loads(open("pass3/dataset.json").read())

## Pick the first root file to get full list of plots
plts = []
f = TFile("pass3/%s/%s.root" % (variation, bkgMCs[0][0]))
for ch in ("el", "mu"):
    moddir = f.Get(ch)
    chdir = moddir.GetDirectory(ch)
    if chdir == None: continue

    for step in [x.GetName() for x in chdir.GetListOfKeys()]:
        stepobj = chdir.Get(step)
        if stepobj == None: continue

        if stepobj.IsA().GetName() in ("TH1D", "TH1F"):
            plts.append({'name':"%s/%s/%s" % (ch, ch, step)})
        elif stepobj.IsA().InheritsFrom("TDirectory"):
            for plt in [x.GetName() for x in stepobj.GetListOfKeys()]:
                if stepobj.Get(plt) == None: continue
                plts.append({'name':"%s/%s/%s/%s" % (ch, ch, step, plt)})

## Start loop
fout = TFile("pass3/%s/preview.root" % variation, "recreate")
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

    dataName = variation+'/'+os.path.basename(fRD[mode].GetName())[:-5]
    #lumi = 1000*sum(x['lumi'] for x in dataset[dataName]['subsamples'])

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
    hMCs = []
    for finName, title, color, f in bkgMCs:
        h = f.Get(plt)
        if h == None: continue
        if 'QCD' in finName: h.Scale(lumi)
        else: h.Scale(lumi)
        hMC.Add(h)
        h.GetStats(stats)
        h.AddBinContent(nbinsX, h.GetBinContent(nbinsX+1))
        h.PutStats(stats)
        h.SetOption("hist")
        h.SetFillColor(color)
        h.SetLineColor(color)
        #h.SetLineStyle(0)
        hsMC.Add(h)
        hMCs.append([title, h])
    hsSig = THStack("hsSig", "")
    for finName, title, color, f in sigMCs:
        h = f.Get(plt)
        h.Scale(lumi)
        h.GetStats(stats)
        h.AddBinContent(nbinsX, h.GetBinContent(nbinsX+1))
        h.PutStats(stats)
        h.SetOption("hist")
        h.SetLineColor(color)
        h.SetLineStyle(1)
        hsSig.Add(h)
    hRatio = hRD.Clone()
    hRatio.Reset()
    hRatio.SetTitle(";%s;Data/MC" % hRD.GetXaxis().GetTitle())
    grpRatio = TGraphErrors()
    rMax = 2
    for b in range(nbinsX):
        yRD, yMC = hRD.GetBinContent(b+1), hMC.GetBinContent(b+1)
        eRD, eMC = hRD.GetBinError(b+1), hMC.GetBinError(b+1)
        r, e = 1e9, 1e9
        if yMC > 0:
            r = yRD/yMC
            rMax = max(r, rMax)
        #if yMC > 0 and yRD > 0: e = r*hypot(eRD/yRD, eMC/yMC)
        if yRD > 0: e = r*abs(eRD/yRD)
        if r == 1e9 or e == 1e9: continue

        x = hRD.GetXaxis().GetBinCenter(b+1)
        w = hRD.GetXaxis().GetBinWidth(b+1)
        n = grpRatio.GetN()
        grpRatio.SetPoint(n, x, r)
        grpRatio.SetPointError(n, w/2, e)
    if rMax > 2: rMax = 3
    hRatio.SetStats(False)
    hRatio.SetMinimum(0.5)
    #hRatio.SetMaximum(rMax)
    hRatio.SetMaximum(1.5)

    ## Draw'em all
    plotDim = (400, 300, 100) # width, main height, ratio height
    margin = (2*40, 20, 60, 60) # left, right, bottom, top
    padH = (plotDim[1] + margin[3], plotDim[2] + margin[2])
    canH = padH[0] + padH[1]
    canW = plotDim[0] + margin[0] + margin[1]

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

    ## Determine allowed maxY for nice drawing
    legY2 = 1.-0.05-1.*margin[3]/padH[0]
    legY1 = legY2-(1+hsMC.GetNhists())*0.05/2
    xmin, xmax = hRD.GetXaxis().GetXmin(), hRD.GetXaxis().GetXmax()
    yMaxOnCanvas = max(hRD.GetMaximum(), hMC.GetMaximum())
    yMaxUnderLeg = 0
    #for ibin in range(hRD.FindBin(xmin+(xmax-xmin)*0.6)-1, hRD.GetNbinsX()): ## ignore the last overflow bin
    for ibin in range(1,hRD.GetNbinsX()): ## ignore the last overflow bin
        yMaxUnderLeg = max(yMaxUnderLeg, hRD.GetBinContent(ibin))
        yMaxUnderLeg = max(yMaxUnderLeg, hMC.GetBinContent(ibin))
    if yMaxUnderLeg/hRD.GetMaximum() < legY1: hRD.SetMaximum(yMaxOnCanvas*1.2)
    else:  hRD.SetMaximum(yMaxOnCanvas*1.2/0.6)

    #pad1.SetLogy()
    #hRD.SetMaximum(hRD.GetMaximum()*100)
    hRD.SetMinimum(0)
    hRD.Draw("")
    hsMC.Draw("same,hist")
    hRD.Draw("same,e")
    hsSig.Draw("same,hist,nostack")
    pad1.RedrawAxis()

    ## Draw legends
    leg = TLegend(0.6, legY1, 1-1.*margin[1]/canW-0.05, legY2, "", "NDC")
    leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    #hists = hsSig.GetHists()
    for title, h in reversed(hMCs):
        leg.AddEntry(h, title, "f")
    leg.AddEntry(hRD, "Data", "lp")
    leg.Draw()

    ## Draw labels
    nLine = 1
    label = TPaveText(1.*margin[0]/canW+0.05, legY2-0.07*nLine, 0.5, legY2, "NDC")
    label.SetBorderSize(0)
    label.SetFillStyle(0)
    label.SetTextAlign(kHAlignLeft+kVAlignBottom)
    chName, stepName = plt.split('/')[0], plt.split('/')[-2]
    label.AddText("%s %s" % (chName, stepName))
    label.Draw()

    fout.cd(dirName)
    pad1.Modified()
    pad1.Update()
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

    for h in (leg, hRD, hMC, hsMC, hRatio, grpRatio, c): del(h)

## Start to print cut flow
cutflow = {
    "count":{"el":{}, "mu":{}},
    "error":{"el":{}, "mu":{}},
    "nstep":0,
    "step":None,
}
nstep = 0
for mode in cutflow["count"].keys():
    h = fRD[mode].Get("%s/%s/cutstep" % (mode, mode))
    nstep = h.GetNbinsX()
    cutflow["count"][mode]["Data"] = [h.GetBinContent(i) for i in range(1, nstep+1)]
    cutflow["error"][mode]["Data"] = [h.GetBinError(i) for i in range(1, nstep+1)]
    if cutflow["step"] == None:
        cutflow["step"] = [h.GetXaxis().GetBinLabel(i) for i in range(1, nstep+1)]

    for finName, title, color, f in bkgMCs:
        h = f.Get("%s/%s/cutstep" % (mode, mode))
        if h == None: continue
        cutflow["count"][mode][finName] = [h.GetBinContent(i) for i in range(1, nstep+1)]
        cutflow["error"][mode][finName] = [h.GetBinError(i) for i in range(1, nstep+1)]

    for finName, title, color, f in sigMCs:
        h = f.Get("%s/%s/cutstep" % (mode, mode))
        if h == None: continue
        cutflow["count"][mode][finName] = [h.GetBinContent(i) for i in range(1, nstep+1)]
        cutflow["error"][mode][finName] = [h.GetBinError(i) for i in range(1, nstep+1)]
cutflow["nstep"] = nstep
printCutflow(cutflow)

## Save cut flow
f = open("pass3/cutflow.json", "w")
f.write(json.dumps(cutflow, indent=4, sort_keys=True))
f.close()

## Save plot list
f = open("pass3/plots.json", "w")
f.write(json.dumps({'plots':plts}, indent=4, sort_keys=True))
f.close()

print "A preview root file for the %s sample is produced" % variation
print "Run `root -l pass3/preview.root' and browse into each directories to open canvases"
print "You can also use dumpRoot command from the hep-tools, dumpRoot pass3/preview.root"
print "Cutflow table is saved in JSON format, pass3/cutflow.json"
