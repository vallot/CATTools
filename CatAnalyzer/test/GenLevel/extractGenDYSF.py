#!/usr/bin/env python

from ROOT import *

xHi = 3*2008.4
xLo20to50 = 3*3205.6-xHi
doAddSF = False
#doAddSF = True

fHi = TFile("DYJets.root")
fLo = TFile("DYJets_10to50.root")
#fHi = TFile("DYJets_MG.root")
#fLo = TFile("DYJets_MG_5to50.root")

hZMHi = fHi.Get("z/default/hZMass")
sumwHi = fHi.Get("z/default/hZPhi").Integral()
hZMHi.Scale(xHi/sumwHi)
hZMHi.SetLineColor(kBlue)
hZMHi.SetFillColor(kBlue)

hZMLo = fLo.Get("z/default/hZMass")
sumwLo = fLo.Get("z/default/hZPhi").Integral()

sumwLo20to50 = hZMLo.Integral(hZMLo.FindBin(20.), hZMLo.FindBin(50-0.001))
sumwLo00to50 = hZMLo.Integral(0, hZMLo.FindBin(50-0.001))
sfLo = sumwLo20to50/sumwLo00to50
hZMLo.Scale(xLo20to50/sumwLo/sfLo)
## Estimate additional factor to patch two region
ym2 = hZMLo.GetBinContent(hZMLo.FindBin(50-0.001)-1)+hZMHi.GetBinContent(hZMHi.FindBin(50-0.001)-1)
ym1 = hZMLo.GetBinContent(hZMLo.FindBin(50-0.001))+hZMHi.GetBinContent(hZMHi.FindBin(50-0.001))
yp1 = hZMLo.GetBinContent(hZMLo.FindBin(50+0.001))+hZMHi.GetBinContent(hZMHi.FindBin(50+0.001))
yp2 = hZMLo.GetBinContent(hZMLo.FindBin(50+0.001)+1)+hZMHi.GetBinContent(hZMHi.FindBin(50+0.001)+1)
rPatch = (ym2/ym1 + yp1/yp2)/2 # Average of ratios between nearby bins
addSF = rPatch/(ym1/yp1)

if doAddSF: hZMLo.Scale(addSF)
hZMLo.SetLineColor(kRed)
hZMLo.SetFillColor(kRed)

hstack = THStack("hstack", "hstack")
hstack.Add(hZMHi)
hstack.Add(hZMLo)

c = TCanvas("c", "c", 500, 500)
c.SetLogy()
hstack.Draw("hist")

print "="*40
print " * Ratio of (Theory)/(Generated) = ", sfLo
print " * Additional factor for patching two region = ", addSF
print "="*40
