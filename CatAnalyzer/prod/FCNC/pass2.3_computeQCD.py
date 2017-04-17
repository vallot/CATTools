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
    "el":[x for x in json.loads(open("pass3/hists.json").read()) if x.startswith("el/el/step")],
    "mu":[x for x in json.loads(open("pass3/hists.json").read()) if x.startswith("mu/mu/step")],
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

## Make mT/dphi histograms to be used in the ABCD method
f = TFile("pass3/nominal/sel_ABCD.root", "recreate")
for ch in ("el", "mu"):
    for dName in dNames[ch]:
        dout = rootGetDir(f, dName)
        for hName in ("kin_mT", "kin_mT_cosDphi"):
            hName = str(dName+'/'+hName)
            hRD     = fRD[ch].Get(hName)
            hRD_QCD = fRD_QCD[ch].Get(hName)
            if None in (hRD, hRD_QCD): continue

            ## Compute (B/D) for each steps
            ##        
            ## QCDSel |  B  |  D
            ## SigSel |  A  |  C 
            ##        +-----+-----
            ##         Iso   NonIso

            dout.cd()
            h_AB = hRD.Clone()
            h_CD = hRD_QCD.Clone()
            h_CD_MC = h_CD.Clone()
            h_AB.SetName(hRD.GetName()+"_AB")
            h_CD.SetName(hRD_QCD.GetName()+"_CD")
            h_CD_MC.Reset()
            h_CD_MC.SetName(hRD_QCD.GetName()+"_CD_MC")
            for title, colour, files in bkgMCs:
                if len(files) < 2: continue
                h_AB.Add(files[0].Get(hName), -lumi)
                h_CD.Add(files[1].Get(hName), -lumi)
                h_CD_MC.Add(files[1].Get(hName), lumi)
            h_AB.Write()
            h_CD.Write()
            h_CD_MC.Write()
f.Close()

"""
## Compute scales, 
scales = {}
for ch in ("el", "mu"):
    for dName in dNames[ch]:
        ## Compute (B/D) for each steps
        ##        
        ## mT>20 |  A  |  C
        ## mT<20 |  B  |  D 
        ##       +-----+-----
        ##         Iso   NonIso

        hName = str(dName+'/kin_mT')
        hRD     = fRD[ch].Get(hName)
        hRD_QCD = fRD_QCD[ch].Get(hName)
        if None in (hRD, hRD_QCD):
            scales[dName] = 1
            continue

        iBinMTCut = hRD.GetXaxis().FindBin(20) ## 20GeV cut
        #nRD_A = hRD.Integral(iBinMTCut, hRD.GetNbinsX()+1) ## High mT, non isolated
        nRD_B = hRD.Integral(0, iBinMTCut-1) ## Low mT, isolated
        #nRD_C = hRD_QCD.Integral(iBinMTCut, hRD.GetNbinsX()+1) ## High mT, non isolated
        nRD_D = hRD_QCD.Integral(0, iBinMTCut-1) ## Low mT, non isolated

        nMC_B, nMC_D = 0., 0.,
        for title, colour, files in bkgMCs:
            if len(files) < 2: continue
            nMC_B += files[0].Get(hName).Integral(0, iBinMTCut-1)
            nMC_D += files[1].Get(hName).Integral(0, iBinMTCut-1)

        scales[dName] = (nRD_B-lumi*nMC_B)/(nRD_D-lumi*nMC_D)

## Do the subtraction to get the shapes
f = TFile("pass3/nominal/QCD_Data.root", "recreate")
for ch in ("el", "mu"):
    for hName in hists[ch]:
        dName = os.path.dirname(hName)
        hName = str(hName)
        scale = scales[dName]

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

"""
