#!/usr/bin/env python

from ROOT import *
import os
from math import sqrt

def combine(fNameCen, fNameUp, fNameDn, fNames, plotNames, combineBy):
    fcen = TFile(fNameCen)
    if fcen == None: return
    fcen.SetBit(TFile.kDevNull)

    print "@@ Processing", fNameCen

    fup = TFile(fNameUp, "recreate")
    fdn = TFile(fNameDn, "recreate")

    files = []
    for fName in fNames:
        files.append(TFile(fName))
        files[-1].SetBit(TFile.kDevNull)

    for pName in plotNames:
        hcen = fcen.Get(pName)
        if hcen == None: continue
        nbins = hcen.GetNbinsX()

        diffs = [[]]*(nbins+2)
        for f in files:
            h = f.Get(pName)
            for b in range(nbins+2):
                diffs[b].append(h.GetBinContent(b)-hcen.GetBinContent(b))

        rdir = os.path.dirname(pName)
        dup = fup.GetDirectory(rdir)
        if dup == None:
            fup.mkdir(rdir)
            dup = fup.GetDirectory(rdir)
        dup.cd()
        hup = hcen.Clone()

        ddn = fdn.GetDirectory(rdir)
        if ddn == None:
            fdn.mkdir(rdir)
            ddn = fdn.GetDirectory(rdir)
        ddn.cd()
        hdn = hcen.Clone()

        if combineBy == "hessian":
            for b in range(nbins+2):
                dqsqr = sum([dyi[b]**2 for dyi in diffs])
                hup.AddBinContent(b, sqrt(dysqr))
                hdn.AddBinContent(b, -sqrt(dysqr))
        elif combineBy == "envelope":
            for b in range(nbins+2):
                dymax = max([dyi[b] for dyi in diffs])
                dymin = min([dyi[b] for dyi in diffs])
                hup.AddBinContent(b, dymax)
                hdn.AddBinContent(b, dymin)
        elif combineBy == "gauss":
            for b in range(nbins+2):
                n = len(diffs)
                dqsqr = sum([dyi[b]**2 for dyi in diffs])
                hup.AddBinContent(b, sqrt(dysqr)/n)
                hdn.AddBinContent(b, -sqrt(dysqr)/n)

        dup.cd()
        hup.Write("", TObject.kOverwrite)

        ddn.cd()
        hdn.Write("", TObject.kOverwrite)

    print "@@ Clean up opened files"
    fcen.Close()
    fup.Close()
    fdn.CLose()
    for f in files:
        gROOT.GetListOfFiles().Remove(f)

    print "@@ Finished", fNameCen

